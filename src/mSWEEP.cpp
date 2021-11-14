#include <iostream>
#include <vector>
#include <exception>
#include <memory>
#include <fstream>

#include "cxxio.hpp"
#include "rcgpar.hpp"

#include "parse_arguments.hpp"
#include "likelihood.hpp"
#include "Sample.hpp"
#include "Reference.hpp"

#include "version.h"
#include "openmp_config.hpp"

#if defined(MSWEEP_OPENMP_SUPPORT) && (MSWEEP_OPENMP_SUPPORT) == 1
#include <omp.h>
#endif

int main (int argc, char *argv[]) {
  std::cerr << "mSWEEP-" << MSWEEP_BUILD_VERSION << " abundance estimation" << std::endl;
  Arguments args;
  if (CmdOptionPresent(argv, argv+argc, "--version")) {
    std::cerr << "mSWEEP-" << MSWEEP_BUILD_VERSION << std::endl;
  }
  if (CmdOptionPresent(argv, argv+argc, "--help")) {
    PrintHelpMessage();
  }
  if (CmdOptionPresent(argv, argv+argc, "--cite")) {
    PrintCitationInfo();
  }
  if (CmdOptionPresent(argv, argv+argc, "--version") || CmdOptionPresent(argv, argv+argc, "--help") || CmdOptionPresent(argv, argv+argc, "--cite")) {
    return 0;
  }
  try {
    ParseArguments(argc, argv, args);
  }
  catch (std::runtime_error &e) {
    std::cerr << "Error in parsing arguments:\n  "
	      << e.what()
	      << "\nexiting" << std::endl;
    return 1;
  }

#if defined(MSWEEP_OPENMP_SUPPORT) && (MSWEEP_OPENMP_SUPPORT) == 1
  omp_set_num_threads(args.optimizer.nr_threads);
#endif

  std::vector<std::unique_ptr<Sample>> samples;
  Reference reference;
  try {
    std::cerr << "Reading the input files" << '\n';
    std::cerr << "  reading group indicators" << '\n';
    if (args.fasta_file.empty()) {
      cxxio::In indicators_file(args.indicators_file);
      reference.read_from_file(indicators_file.stream(), args.groups_list_delimiter);
    } else {
      cxxio::In groups_file(args.groups_list_file);
      cxxio::In fasta_file(args.fasta_file);
      reference.match_with_fasta(args.groups_list_delimiter, groups_file.stream(), fasta_file.stream());
    }

    std::cerr << "  read " << reference.get_n_refs() << " group indicators" << std::endl;

    std::cerr << (args.read_likelihood_mode ? "  reading likelihoods from file" : "  reading pseudoalignments") << '\n';
    if (args.bootstrap_mode) {
      samples.emplace_back(new BootstrapSample(args.seed));
    } else {
      samples.emplace_back(new Sample());
    }
    if (!args.themisto_mode && !args.read_likelihood_mode) {
      // Check that the number of reference sequences matches in the grouping and the alignment.
      reference.verify_kallisto_alignment(*args.infiles.run_info);
      ReadKallisto(reference.get_n_refs(), *args.infiles.ec, *args.infiles.tsv, &samples.back()->pseudos);
    } else if (!args.read_likelihood_mode) {
      if (!args.themisto_index_path.empty()) {
	cxxio::In themisto_index(args.themisto_index_path + "/coloring-names.txt");
	reference.verify_themisto_index(themisto_index);
      }
      cxxio::In forward_strand(args.tinfile1);
      cxxio::In reverse_strand(args.tinfile2);
      std::vector<std::istream*> strands = { &forward_strand.stream(), &reverse_strand.stream() };
      ReadThemisto(get_mode(args.themisto_merge_mode), reference.get_n_refs(), strands, &samples.back()->pseudos);

    } else {
      if (reference.get_n_groupings() > 1) {
	throw std::runtime_error("Using more than one grouping with --read-likelihood is not yet implemented.");
      }
      cxxio::In likelihoods(args.likelihood_file);
      samples.back()->read_likelihood(reference.get_grouping(0), likelihoods.stream());
    }
    samples.back()->process_aln(args.bootstrap_mode);

    std::cerr << "  read " << (args.batch_mode ? samples.size() : samples[0]->num_ecs()) << (args.batch_mode ? " samples from the batch" : " unique alignments") << std::endl;
  } catch (std::runtime_error &e) {
    std::cerr << "Reading the input files failed:\n  ";
    std::cerr << e.what();
    std::cerr << "\nexiting" << std::endl;
    return 1;
  }


  // Estimate abundances with all groupings that were provided
  uint16_t n_groupings = reference.get_n_groupings();
  std::string outfile_name = args.outfile;
  for (uint16_t i = 0; i < n_groupings; ++i) {
    const Grouping *grouping = &reference.get_grouping(i);
    uint32_t n_groups = grouping->get_n_groups();

    // Initialize the prior counts on the groups
    args.optimizer.alphas = std::vector<double>(n_groups, 1.0);

    args.outfile = outfile_name;
    // Set output file name correctly
    if (n_groupings > 1 && !outfile_name.empty()) { // Backwards compatibility with v1.4.0 or older in output names
      args.outfile += "_";
      args.outfile += std::to_string(i);
    }

    std::cerr << "Building log-likelihood array" << std::endl;
    for (uint16_t j = 0; j < samples.size(); ++j) {
      if (!args.read_likelihood_mode) {
	likelihood_array_mat((*grouping), reference.get_group_indicators(i), args.optimizer.bb_constants, (*samples[j]));
	if (j == n_groupings - 1) {
	  // Free memory used by the configs after all likelihood matrices are built.
	  samples[j]->pseudos.ec_configs.clear();
	  samples[j]->pseudos.ec_configs.shrink_to_fit();
	}
      }

      if (args.optimizer.write_likelihood || args.optimizer.write_likelihood_bitseq) {
	cxxio::Out of;
	std::string outfile(args.outfile);
	outfile += (args.optimizer.write_likelihood_bitseq ? "_bitseq" : "");
	outfile += "likelihoods.txt";
	if (args.optimizer.gzip_probs) {
	  outfile += ".gz";
	  of.open_compressed(outfile);
	} else {
	  of.open(outfile);
	}
	if (args.optimizer.write_likelihood_bitseq) {
	  std::cerr << "Writing likelihood matrix in BitSeq format" << std::endl;
	  samples[j]->write_likelihood_bitseq(args.optimizer.gzip_probs, n_groups, of.stream());
	} else {
	  std::cerr << "Writing likelihood matrix" << std::endl;
	  samples[j]->write_likelihood(args.optimizer.gzip_probs, n_groups, of.stream());
	}
      }
    }

    // Process the reads accordingly
    if (args.optimizer.no_fit_model) {
      std::cerr << "Skipping relative abundance estimation (--no-fit-model toggled)" << std::endl;
    } else {
      std::cerr << "Estimating relative abundances" << std::endl;
      for (uint32_t i = 0; i < samples.size(); ++i) {
	// Set out buffer
	std::string outfile(args.outfile);
	outfile = (outfile.empty() || !args.batch_mode ? outfile : outfile + '/' + samples[i]->cell_name());
	cxxio::Out of;
	if (!outfile.empty()) {
	  outfile += "_abundances.txt";
	  of.open(outfile);
	}

	if (args.bootstrap_mode) {
	  std::cerr << "Running estimation with " << args.iters << " bootstrap iterations" << '\n';
	  BootstrapSample* bs = static_cast<BootstrapSample*>(&(*samples[i]));
	  bs->bootstrap_abundances((*grouping), args);
	  bs->write_bootstrap(grouping->get_names(), args.iters, args.batch_mode, (outfile.empty() ? std::cout : of.stream()));
	} else {
	  // Estimate relative abundances
	  samples[i]->ec_probs = rcgpar::rcg_optl_omp(samples[i]->ll_mat, samples[i]->log_ec_counts, args.optimizer.alphas, args.optimizer.tolerance, args.optimizer.max_iters, std::cerr);
	  samples[i]->relative_abundances = rcgpar::mixture_components(samples[i]->ec_probs, samples[i]->log_ec_counts);
	  samples[i]->write_abundances(grouping->get_names(), (outfile.empty() ? std::cout : of.stream()));
	}
	// Write the probability matrix
	std::string probs_outfile(args.outfile);
	if ((args.optimizer.write_probs || args.optimizer.print_probs) && !args.outfile.empty()) {
	  cxxio::Out of;
	  probs_outfile += "_probs.csv";
	  if (args.optimizer.gzip_probs) {
	    probs_outfile += ".gz";
	    of.open_compressed(probs_outfile);
	  } else if (!args.optimizer.print_probs){
	    of.open(probs_outfile);
	  }
	  samples[i]->write_probabilities(grouping->get_names(), (args.optimizer.print_probs ? std::cout : of.stream()));
	}
      }
    }
  }
  return 0;
}

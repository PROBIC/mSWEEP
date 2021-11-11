#include <iostream>
#include <vector>
#include <exception>
#include <memory>

#include "bxzstr.hpp"
#include "cxxio.hpp"

#include "parse_arguments.hpp"
#include "read_pseudoalignment.hpp"
#include "process_reads.hpp"
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
    if (!args.themisto_mode && !args.read_likelihood_mode) {
      // Check that the number of reference sequences matches in the grouping and the alignment.
      reference.verify_kallisto_alignment(*args.infiles.run_info);
      ReadPseudoalignment(args.infiles, reference.get_n_refs(), samples, args.bootstrap_mode);
    } else if (!args.read_likelihood_mode) {
      if (!args.themisto_index_path.empty()) {
	cxxio::In themisto_index(args.themisto_index_path + "/coloring-names.txt");
	reference.verify_themisto_index(themisto_index);
      }
      ReadPseudoalignment(args.tinfile1, args.tinfile2, args.themisto_merge_mode, args.bootstrap_mode, reference.get_n_refs(), samples);
    } else {
      if (reference.get_n_groupings() > 1) {
	throw std::runtime_error("Using more than one grouping with --read-likelihood is not yet implemented.");
      }
      cxxio::In likelihoods(args.likelihood_file);
      ReadLikelihood(reference.get_grouping(0), args.bootstrap_mode, likelihoods.stream(), samples);
    }

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
    uint32_t n_groups = reference.get_grouping(i).get_n_groups();

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
	samples[j]->CalcLikelihood(reference.get_grouping(i), args.optimizer.bb_constants, reference.get_group_indicators(i), n_groupings == 1);
      }

      if (args.optimizer.write_likelihood || args.optimizer.write_likelihood_bitseq) {
	std::cerr << "Writing likelihood matrix" << std::endl;
	if (args.optimizer.write_likelihood) {
	  samples[j]->write_likelihood(args.optimizer.gzip_probs, n_groups, args.outfile);
	}
	if (args.optimizer.write_likelihood_bitseq) {
	  samples[j]->write_likelihood_bitseq(args.optimizer.gzip_probs, n_groups, args.outfile);
	}
      }
    }

    // Process the reads accordingly
    if (args.optimizer.no_fit_model) {
      std::cerr << "Skipping relative abundance estimation (--no-fit-model toggled)" << std::endl;
    } else {
      std::cerr << "Estimating relative abundances" << std::endl;
      if (args.bootstrap_mode) {
	ProcessBootstrap(reference.get_grouping(i), args, samples);
      } else {
	ProcessReads(reference.get_grouping(i), args, samples);
      }
    }
  }

  return 0;
}

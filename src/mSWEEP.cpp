#include "mSWEEP.hpp"

#include <iostream>
#include <exception>
#include <fstream>

#include "cxxio.hpp"
#include "rcgpar.hpp"

#include "likelihood.hpp"
#include "version.h"

void ReadInput(const Arguments &args, std::vector<std::unique_ptr<Sample>> *samples, std::ostream &log, Reference *reference) {
  log << "Reading the input files" << '\n';
  log << "  reading group indicators" << '\n';
  if (args.fasta_file.empty()) {
    cxxio::In indicators_file(args.indicators_file);
    reference->read_from_file(indicators_file.stream(), args.groups_list_delimiter);
  } else {
    cxxio::In groups_file(args.groups_list_file);
    cxxio::In fasta_file(args.fasta_file);
    reference->match_with_fasta(args.groups_list_delimiter, groups_file.stream(), fasta_file.stream());
  }

  log << "  read " << reference->get_n_refs() << " group indicators" << '\n';

  log << (args.read_likelihood_mode ? "  reading likelihoods from file" : "  reading pseudoalignments") << '\n';
  if (!args.themisto_mode && !args.read_likelihood_mode) {
    // Check that the number of reference sequences matches in the grouping and the alignment.
    reference->verify_kallisto_alignment(*args.infiles.run_info);
    ReadKallisto(reference->get_n_refs(), *args.infiles.ec, *args.infiles.tsv, &samples->back()->pseudos);
  } else if (!args.read_likelihood_mode) {
    if (!args.themisto_index_path.empty()) {
      cxxio::In themisto_index(args.themisto_index_path + "/coloring-names.txt");
      reference->verify_themisto_index(themisto_index);
    }
    cxxio::In forward_strand(args.tinfile1);
    cxxio::In reverse_strand(args.tinfile2);
    std::vector<std::istream*> strands = { &forward_strand.stream(), &reverse_strand.stream() };
    ReadThemisto(get_mode(args.themisto_merge_mode), reference->get_n_refs(), strands, &samples->back()->pseudos);
  } else {
    if (reference->get_n_groupings() > 1) {
      throw std::runtime_error("Using more than one grouping with --read-likelihood is not yet implemented.");
    }
    cxxio::In likelihoods(args.likelihood_file);
    samples->back()->read_likelihood(reference->get_grouping(0), likelihoods.stream());
  }
  samples->back()->process_aln(args.bootstrap_mode);

  log << "  read " << (args.batch_mode ? samples->size() : (*samples)[0]->num_ecs()) << (args.batch_mode ? " samples from the batch" : " unique alignments") << '\n';
  log.flush();
}

void ConstructLikelihood(const Arguments &args, const Grouping &grouping, const std::vector<uint32_t> &group_indicators, const std::unique_ptr<Sample> &sample, bool free_ec_counts) {
  if (!args.read_likelihood_mode) {
    likelihood_array_mat(grouping, group_indicators, args.optimizer.bb_constants, (*sample));
  }
  if (free_ec_counts) {
    // Free memory used by the configs after all likelihood matrices are built.
    sample->pseudos.ec_configs.clear();
    sample->pseudos.ec_configs.shrink_to_fit();
  }
}

void WriteResults(const Arguments &args, const std::unique_ptr<Sample> &sample, const Grouping &grouping, const uint16_t n_groupings, const uint16_t current_grouping) {
  cxxio::Out of;

  // Set output file name correctly
  // for backwards compatibility with v1.4.0 or older
  std::string outfile = args.outfile;
  if (n_groupings > 1 && !args.outfile.empty()) {
    outfile += "_";
    outfile += std::to_string(current_grouping);
  }

  // Write likelihoods
  if (args.optimizer.write_likelihood || args.optimizer.write_likelihood_bitseq) {
    std::string ll_outfile(outfile);
    ll_outfile += (args.optimizer.write_likelihood_bitseq ? "_bitseq" : "");
    ll_outfile += "likelihoods.txt";
    if (args.optimizer.gzip_probs) {
      ll_outfile += ".gz";
      of.open_compressed(ll_outfile);
    } else {
      of.open(ll_outfile);
    }
    if (args.optimizer.write_likelihood_bitseq) {
      sample->write_likelihood_bitseq(grouping.get_n_groups(), of.stream());
    } else {
      sample->write_likelihood(grouping.get_n_groups(), of.stream());
    }
  }

  // Relative abundances
  std::string abundances_outfile(outfile);
  abundances_outfile = (args.outfile.empty() || !args.batch_mode ? abundances_outfile : abundances_outfile + '/' + sample->cell_name());
  if (!args.outfile.empty()) {
    abundances_outfile += "_abundances.txt";
    of.open(abundances_outfile);
  }
  if (args.bootstrap_mode) {
    BootstrapSample* bs = static_cast<BootstrapSample*>(&(*sample));
    bs->write_bootstrap(grouping.get_names(), args.iters, (args.outfile.empty() ? std::cout : of.stream()));
  } else {
    sample->write_abundances(grouping.get_names(), (args.outfile.empty() ? std::cout : of.stream()));
  }

  // Probability matrix
  std::string probs_outfile(outfile);
  if ((args.optimizer.write_probs || args.optimizer.print_probs) && !args.outfile.empty()) {
    probs_outfile += "_probs.csv";
    if (args.optimizer.gzip_probs) {
      probs_outfile += ".gz";
      of.open_compressed(probs_outfile);
    } else if (!args.optimizer.print_probs){
      of.open(probs_outfile);
    }
    sample->write_probabilities(grouping.get_names(), (args.optimizer.print_probs ? std::cout : of.stream()));
  }
}

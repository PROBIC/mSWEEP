#include "mSWEEP.hpp"

#include <iostream>
#include <exception>
#include <fstream>

#include "cxxio.hpp"
#include "rcgpar.hpp"

#include "likelihood.hpp"
#include "version.h"

void VerifyThemistoIndex(const std::string &themisto_index_path, const Reference &reference) {
  // TODO merge with reference.verify_themisto_index
  if (!themisto_index_path.empty()) {
    try {
      cxxio::In themisto_index(themisto_index_path + "/coloring-names.txt");
      reference.verify_themisto_index(themisto_index);
    } catch (const std::runtime_error &e) {
      throw std::runtime_error("--themisto-index flag is not supported for Themisto v2.0.0 or newer:\n" + std::string(e.what()));
    } catch (const std::domain_error &e) {
      throw e;
    }
  }
}

void ReadGroupIndicators(const Arguments &args, Reference *reference) {
  if (args.fasta_file.empty()) {
    cxxio::In indicators_file(args.indicators_file);
    reference->read_from_file(indicators_file.stream(), args.groups_list_delimiter);
  } else {
    cxxio::In groups_file(args.groups_list_file);
    cxxio::In fasta_file(args.fasta_file);
    reference->match_with_fasta(args.groups_list_delimiter, groups_file.stream(), fasta_file.stream());
  }
}

void ReadPseudoalignments(const Arguments &args, const Reference &reference, std::unique_ptr<Sample> &sample) {
  sample.reset(new Sample(reference)); // For some reason the sample needs to be reset here ??
  VerifyThemistoIndex(args.themisto_index_path, reference);
  cxxio::In forward_strand(args.tinfile1);
  cxxio::In reverse_strand(args.tinfile2);
  std::vector<std::istream*> strands = { &forward_strand.stream(), &reverse_strand.stream() };

  if (args.compact_alignments) {
    sample->pseudos.set_parse_from_buffered();
  }
  telescope::read::ThemistoGrouped(telescope::get_mode(args.themisto_merge_mode), strands, &sample->pseudos);
}

void ReadLikelihoodFromFile(const Arguments &args, const Reference &reference, std::ostream &log, const std::unique_ptr<Sample> &sample) {
  log << "  reading likelihoods from file" << '\n';
  if (reference.get_n_groupings() > 1) {
    throw std::runtime_error("Using more than one grouping with --read-likelihood is not yet implemented.");
  }
  cxxio::In likelihoods(args.likelihood_file);
  sample->read_likelihood(reference.get_grouping(0), likelihoods.stream());
}

void ConstructLikelihood(const Arguments &args, const Grouping &grouping, const std::vector<uint32_t> &group_indicators, const std::unique_ptr<Sample> &sample, bool free_ec_counts) {
  if (!args.read_likelihood_mode) {
    likelihood_array_mat(grouping, group_indicators, args.optimizer.bb_constants, (*sample));
  }
  if (free_ec_counts) {
    // Free memory used by the configs after all likelihood matrices are built.
    //    sample->pseudos.clear_configs();
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
    ll_outfile += "_likelihoods.txt";
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
  if (!args.optimizer.no_fit_model) {
    std::string abundances_outfile(outfile);
    if (!args.outfile.empty()) {
      std::string abundances_outfile = outfile + "_abundances.txt";
      of.open(abundances_outfile);
    }
    if (args.bootstrap_mode) {
      BootstrapSample* bs = static_cast<BootstrapSample*>(&(*sample));
      bs->write_bootstrap(grouping.get_names(), args.iters, (args.outfile.empty() ? std::cout : of.stream()));
    } else {
      sample->write_abundances(grouping.get_names(), (args.outfile.empty() ? std::cout : of.stream()));
    }
  }

  // Probability matrix
  if (args.optimizer.print_probs && !args.optimizer.no_fit_model) {
    sample->write_probabilities(grouping.get_names(), std::cout);
  }
  if (args.optimizer.write_probs && !args.optimizer.no_fit_model) {
    std::string probs_outfile(outfile);
    probs_outfile += "_probs.csv";
    if (args.optimizer.gzip_probs) {
      probs_outfile += ".gz";
      of.open_compressed(probs_outfile);
    } else {
      of.open(probs_outfile);
    }
    sample->write_probabilities(grouping.get_names(), (outfile.empty() ? std::cout : of.stream()));
  }
}

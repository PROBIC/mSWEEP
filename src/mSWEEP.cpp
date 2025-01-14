// mSWEEP: Estimate abundances of reference lineages in DNA sequencing reads.
//
// MIT License
//
// Copyright (c) 2023 Probabilistic Inference and Computational Biology group @ UH
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <cstddef>
#include <fstream>
#include <exception>
#include <memory>
#include <limits>

#include "cxxargs.hpp"
#include "Matrix.hpp"
#include "cxxio.hpp"
#include "rcgpar.hpp"
#include "bin_reads.h"

#include "mSWEEP_alignment.hpp"

#include "mSWEEP_openmp_config.hpp"
#include "mSWEEP_version.h"

#include "mSWEEP_log.hpp"
#include "Reference.hpp"
#include "Sample.hpp"
#include "Likelihood.hpp"
#include "OutfileDesignator.hpp"

void PrintCitationInfo() {
  // Print the citation information to cerr.
  // TODO: add information about citing mGEMS when binning.
  std::cerr << "Please cite us as:\n"
	    << "\tMäklin T, Kallonen T, David S et al. High-resolution sweep\n"
	    << "\tmetagenomics using fast probabilistic inference [version 2;\n"
	    << "\tpeer review: 2 approved]. Wellcome Open Res 2021, 5:14\n"
	    << "\t(https://doi.org/10.12688/wellcomeopenres.15639.2)" << std::endl;
}

bool CmdOptionPresent(char **begin, char **end, const std::string &option) {
  // Check if `option` is present the char array pointed to by `begin` and `end`.
  // Returns `true` or `false` for is present / is not present.
  return (std::find(begin, end, option) != end);
}

void parse_args(int argc, char* argv[], cxxargs::Arguments &args) {
  // Parse the command line arguments using the cxxargs library.
  args.add_long_argument<bool>("verbose", "Print status messages to cerr.", false);
  args.add_long_argument<bool>("version", "Print mSWEEP version.", false);
  args.add_long_argument<bool>("cite", "Print citation information.", false);
  args.add_long_argument<bool>("help", "Print the help message.\n\nPseudoalignment files (required: -1 and -2, or only -x; will try to read from cin if none are given):", false);

  // Pseudoalignment files
  args.add_long_argument<std::string>("themisto-1", "Pseudoalignment results from Themisto for the 1st strand of paired-end reads.");
  args.add_long_argument<std::string>("themisto-2", "Pseudoalignment results from Themisto for the 2nd strand of paired-end reads.");

  // Separate pseudoalignment files are not required if supplied via a list
  args.add_long_argument<std::vector<std::string>>("themisto", "Single themisto alignment file or a comma separated list of several files.\n\nGroup indicators (required):");
  args.set_not_required("themisto");
  args.set_not_required("themisto-1");
  args.set_not_required("themisto-2");

  // Cluster indicators
  args.add_short_argument<std::string>('i', "Group indicators for the pseudoalignment reference.\n\nOutput prefix:");

  // Output prefix
  args.add_short_argument<std::string>('o', "Prefix for output files written from mSWEEP (default: print to cout).\n\nBinning options:", "");

  // Run the mGEMS binning algorithm
  args.add_long_argument<bool>("bin-reads", "Run the mGEMS binning algorithm and write bins to the directory `-o` points to (default: false).", false);
  args.add_long_argument<std::vector<std::string>>("target-groups", "Only extract these groups, supply as comma separated list (default: extract all groups).");
  args.set_not_required("target-groups");
  args.set_not_required("target-groups");
  args.add_long_argument<double>("min-abundance", "Only extract groups that have a relative abundance higher than this value (default: 0).\n\nOutput options:");
  args.set_not_required("min-abundance");

  // Options for outputting the probability matrix
  args.add_long_argument<bool>("write-probs", "If specified, write the estimated read-to-group probabilities to a file with \"_probs.tsv\" suffix (default:false).", false);
  args.add_long_argument<bool>("print-probs", "Print the read equivalence class probabilities to cout even if `-o` is given (default: false).", false);

  // Write likelihood or not
  args.add_long_argument<bool>("write-likelihood", "Write the internal likelihood matrix to a file with \"_likelihoods.txt\" suffix (default: false).", false);
  args.add_long_argument<bool>("write-likelihood-bitseq", "Write the likelihoods in a format can be parsed by BitSeq's (https://github.com/bitseq/bitseq) functions (default: false).", false);

  // Toggle compression
  args.add_long_argument<std::string>("compress", "Compress all output files using the given algorithm (one of z, bz2, lzma; default: don't compress).", "plaintext");
  args.add_long_argument<int>("compression-level", "Compression level (0-9; default: 6).\n\nInput options:", 6);

  // How to merge paired alignments
  args.add_long_argument<std::string>("themisto-mode", "How to merge pseudoalignments for paired-end reads (intersection, union, or unpaired; default: intersection).", "intersection");

  // Pseudoalignments are not required if reading likelihood from a file
  args.add_long_argument<std::string>("read-likelihood", "Path to a precomputed likelihood file written with the --write-likelihood toggle. Can't be used with --bin-reads.\n\nEstimation options:");
  args.set_not_required("read-likelihood");

  // Number of threads for parallel estimation
  args.add_short_argument<size_t>('t', "How many threads to use in abundance estimation (default: 1).", (size_t)1);
  // Skip running the optimizer (useful for just writing the likelihood matrix)
  args.add_long_argument<bool>("no-fit-model", "Do not estimate the abundances. Useful if only the likelihood matrix is required (default: false).", false);
  // Maximum iterations to run the optimizer for
  args.add_long_argument<size_t>("max-iters", "Maximum number of iterations to run the abundance estimation optimizer for (default: 5000).", (size_t)5000);
  // Tolerance for abundance estimation convergence
  args.add_long_argument<double>("tol", "Optimization terminates when the bound changes by less than the given tolerance (default: 0.000001).", (double)0.000001);
  // Algorithm to use for abundance estimation
  args.add_long_argument<std::string>("algorithm", "Which algorithm to use for abundance estimation (one of rcggpu, emgpu, rcgcpu (original mSWEEP); default: rcgcpu).", "rcgcpu");
  // Precision for abundance estimation with emgpu algorithm
  args.add_long_argument<std::string>("emprecision", "Precision to use for the emgpu algorithm (one of float, double; default: double).\n\nBootstrapping options:", "double");

  // Number of iterations to run bootstrapping for
  args.add_long_argument<size_t>("iters", "Number of times to rerun estimation with bootstrapped alignments (default: 0).", (size_t)0);
  // Seed for bootstrapping
  args.add_long_argument<size_t>("seed", "Seed for the random generator used in bootstrapping (default: random).", 26012023);
  // How many reads to resample when bootstrapping
  args.add_long_argument<size_t>("bootstrap-count", "How many pseudoalignments to resample when bootstrapping (default: number of reads).\n\nLikelihood options:", 0);

  // Mean fraction of aligned sequences for the likelihood
  args.add_short_argument<double>('q', "Mean for the beta-binomial component (default: 0.65).", 0.65);
  // Dispersion term for likelihood
  args.add_short_argument<double>('e', "Dispersion term for the beta-binomial component (default: 0.01).", 0.01);
  // Prior parameters for estimation
  args.add_long_argument<std::vector<double>>("alphas", "Prior counts for the relative abundances, supply as comma-separated nonzero values (default: all 1.0).");
  args.add_long_argument<double>("zero-inflation", "Likelihood of an observation that contains 0 pseudoalignments against a reference group (default: 0.01).\n\nExperimental options:", 0.01);
  args.set_not_required("alphas");

  args.add_long_argument<bool>("run-rate", "Calculate relative reliability for each abundance estimate using RATE (default: false).", false);
  args.add_long_argument<size_t>("min-hits", "Only consider target groups that have at least this many reads align to any sequence in them (default: 0).", (size_t)0);

  if (CmdOptionPresent(argv, argv+argc, "--help")) {
    // Print help message and continue.
    std::cerr << "\n" + args.help() << '\n' << '\n';
  }
  args.parse(argc, argv);

  if (!CmdOptionPresent(argv, argv+argc, "--themisto") && CmdOptionPresent(argv, argv+argc, "--themisto-1") && CmdOptionPresent(argv, argv+argc, "--themisto-2")) {
    // Value of "themisto" is used internally to access the reads.
    args.set_val<std::vector<std::string>>("themisto", std::vector<std::string>({ args.value<std::string>("themisto-1"), args.value<std::string>("themisto-2") }));
  }
}

void finalize(const std::string &msg, mSWEEP::Log &log, bool abort = false) {
  // Set the state of the program so that it can finish correctly:
  // - Writes a potential message `msg` to the log.
  // - Flushes the log (ensure all messages are displayed).
  //
  // Input:
  //   `msg`: message to print.
  //   `log`: logger (see msweep_log.hpp).
  //
  if (!msg.empty())
    std::cerr << msg;
  log.flush();
}

seamat::DenseMatrix<double> rcg_optl(const cxxargs::Arguments &args, const seamat::Matrix<double> &ll_mat, const std::vector<double> &log_ec_counts, const std::vector<double> &prior_counts, mSWEEP::Log &log) {
  // Wrapper for calling rcgpar with omp or gpu depending on config.
  //
  // Input:
  //   `args`: commandl line arguments.
  //   `ll_mat`: the likelihood matrix constructed with functions in likelihood.cpp.
  //   `log_ec_counts`: natural logarithms of read equivalence class observation counts.
  //   `prior_counts`: prior counts for the mixture model components.
  //   `log`: logger (see msweep_log.hpp).
  //
  // Output:
  //   `ec_probs`: read-lineage probability matrix. Use rcgpar::mixture_components to
  //               transform this into the relative abundances vector.
  //
  std::ofstream of;

  if (args.value<std::string>("algorithm") == "rcggpu") {
    // Run rcg on CPU or GPU if present
    const seamat::DenseMatrix<double> &ec_probs = rcgpar::rcg_optl_torch(ll_mat, log_ec_counts, prior_counts, args.value<double>("tol"), args.value<size_t>("max-iters"), (args.value<bool>("verbose") ? log.stream() : of));
    return ec_probs;
  } else if (args.value<std::string>("algorithm") == "rcgcpu") {
    // Run rcg on CPU, use OpenMP if supported
    const seamat::DenseMatrix<double> &ec_probs = rcgpar::rcg_optl_omp(ll_mat, log_ec_counts, prior_counts, args.value<double>("tol"), args.value<size_t>("max-iters"), (args.value<bool>("verbose") ? log.stream() : of));
    return ec_probs;
  } else {
    // Run em on CPU or GPU if present
    const seamat::DenseMatrix<double> &ec_probs = rcgpar::em_torch(ll_mat, log_ec_counts, prior_counts, args.value<double>("tol"), args.value<size_t>("max-iters"), (args.value<bool>("verbose") ? log.stream() : of), args.value<std::string>("emprecision"));
    return ec_probs;
  }
}

int main (int argc, char *argv[]) {
  // mSWEEP executable main
  
  mSWEEP::Log log(std::cerr, CmdOptionPresent(argv, argv+argc, "--verbose"), false); // logger class from msweep_log.hpp

  // Use the cxxargs library to parse the command line arguments.
  cxxargs::Arguments args("mSWEEP", "Usage: mSWEEP --themisto-1 <forwardPseudoalignments> --themisto-2 <reversePseudoalignments> -i <groupIndicatorsFile>");

  // Print version when running if `--verbose` is used.
  log << "mSWEEP-" << MSWEEP_BUILD_VERSION << " abundance estimation" << '\n';
  if (CmdOptionPresent(argv, argv+argc, "--version")) {
    // Print the version if explicitly requested (ignores `--verbose`).
    log.status(std::string("mSWEEP-") + std::string(MSWEEP_BUILD_VERSION));
  }
  if (CmdOptionPresent(argv, argv+argc, "--cite")) {
    // Print the citation information.
    PrintCitationInfo();
  }

  // Actually parse the arguments here
  try {
    // cxxargs throws an error if required arguments are not supplied, catch below.
    log << "Parsing arguments" << '\n';
    parse_args(argc, argv, args);

    // If a output file was supplied check that the directory it's in exists
    if (args.value<std::string>('o').find("/") != std::string::npos) {
      // Assume that if printing to the current directory it exists.
      std::string outfile_dir = args.value<std::string>('o');
      outfile_dir.erase(outfile_dir.rfind("/"), outfile_dir.size());
      cxxio::directory_exists(outfile_dir);
    }

  } catch (const std::exception &e) {
    // TODO: catch the different cxxargs exception types and print informative error messages.
    finalize("Error in parsing arguments:\n  " + std::string(e.what()) + "\nexiting\n", log);
    finalize("", log);
    return 1;
  }

  // Exit if any of `--help`, `--version` or `--cite` were given.
  if (CmdOptionPresent(argv, argv+argc, "--version") || CmdOptionPresent(argv, argv+argc, "--help") || CmdOptionPresent(argv, argv+argc, "--cite")) {
    finalize("", log);
    return 0;
  }

  // Get the number of threads from the arguments and configure OpenMP if enabled.
#if defined(MSWEEP_OPENMP_SUPPORT) && (MSWEEP_OPENMP_SUPPORT) == 1
  omp_set_num_threads(args.value<size_t>('t'));
#endif

  // First read the group indicators
  std::unique_ptr<mSWEEP::Reference> reference;
  log << "Reading the input files" << '\n';
  try {
      log << "  reading group indicators" << '\n';
      cxxio::In indicators(args.value<std::string>('i'));
      reference = std::move(mSWEEP::ConstructAdaptiveReference(&indicators.stream(), '\t'));
      if (reference->get_n_groupings() > 1) {
	log << "  read " << reference->get_n_groupings() << " groupings" << '\n';
      }
      log << "  read " << reference->get_n_refs() << " group indicators" << '\n';
  } catch (std::exception &e) {
    finalize("Reading group indicators failed:\n  " + std::string(e.what()) + "\nexiting\n", log, true);
    return 1;
  }

  // Estimate abundances with all groupings that were provided
  uint16_t n_groupings;
  n_groupings = reference->get_n_groupings();

  // Wrapper class for ensuring the outfile names are set consistently and correctly
  mSWEEP::OutfileDesignator out(args.value<std::string>('o'), n_groupings, args.value<std::string>("compress"), args.value<int>("compression-level"));

  // Estimate abundances with all groupings
  for (uint16_t i = 0; i < n_groupings; ++i) {
      size_t n_groups = reference->n_groups(i);

      // `Sample` and its children are classes for storing data from
      // the pseudoalignment that are needed or not needed depending on
      // the command line arguments.
      std::unique_ptr<mSWEEP::Sample> sample = std::make_unique<mSWEEP::PlainSample>();
      bool bootstrap_mode = args.value<size_t>("iters") > (size_t)0;
      bool bin_reads = args.value<bool>("bin-reads");

      // These are the main inputs to the abundance estimation code.
      size_t n_refs = reference->get_n_refs();
      std::unique_ptr<mSWEEP::Likelihood<double>> log_likelihoods;

      // Check if reading likelihood from file.
      if (!CmdOptionPresent(argv, argv+argc, "--read-likelihood")) {
	// Start from the pseudoalignments.
	// To save memory, the alignment can go out of scope.
	// The necessary values are stored in the Sample class.
	log << "  reading pseudoalignments" << '\n';
	mSWEEP::Alignment alignment(n_refs);
	// std::unique_ptr<telescope::Alignment> alignment;
	try {
	  const std::vector<std::string> &alignment_paths = args.value<std::vector<std::string>>("themisto");
	  size_t n_files = alignment_paths.size();

	  // Open the infiles
	  std::vector<cxxio::In> infiles;
	  infiles.reserve(n_files);
	  std::vector<std::istream*> strands(n_files);
	  if (n_files > 0) {
	    for (size_t i = 0; i < n_files; ++i) {
	      infiles.emplace_back(cxxio::In(alignment_paths[i]));
	      strands[i] = &infiles[i].stream();
	    }
	  } else {
	    strands.emplace_back(&std::cin);
	  }

	  alignment.read(args.value<std::string>("themisto-mode"), strands);
	  log << "  read alignments for " << alignment.n_reads() << " reads" << '\n';
	  log << "Building equivalence classes" << '\n';
	  alignment.collapse();
	  log << "  found " << alignment.n_ecs() << " unique alignments" << '\n';

	} catch (std::exception &e) {
	  finalize("Reading the pseudoalignments failed:\n  " + std::string(e.what()) + "\nexiting\n", log, true);
	  return 1;
	}

	// Use the alignment data to populate the log_likelihoods matrix.
	try {
	    // Read the pseudoalignment
	    log << "Computing the likelihood matrix" << '\n';
	    if (n_groups <= std::numeric_limits<uint8_t>::max()) {
		alignment.add_groups(static_cast<const mSWEEP::AdaptiveReference<uint8_t>*>(&(*reference))->get_group_indicators(i));
	    } else if (n_groups <= std::numeric_limits<uint16_t>::max()) {
		alignment.add_groups(static_cast<const mSWEEP::AdaptiveReference<uint16_t>*>(&(*reference))->get_group_indicators(i));
	    } else if (n_groups <= std::numeric_limits<uint32_t>::max()) {
		alignment.add_groups(static_cast<const mSWEEP::AdaptiveReference<uint32_t>*>(&(*reference))->get_group_indicators(i));
	    } else {
		alignment.add_groups(static_cast<const mSWEEP::AdaptiveReference<uint64_t>*>(&(*reference))->get_group_indicators(i));
	    }

	    log_likelihoods = mSWEEP::ConstructAdaptiveLikelihood<double>(alignment, reference->get_grouping(i), args.value<double>('q'), args.value<double>('e'), args.value<size_t>("min-hits"), args.value<double>("zero-inflation"));
	}  catch (std::exception &e) {
	  finalize("Building the log-likelihood array failed:\n  " + std::string(e.what()) + "\nexiting\n", log, true);
	  return 1;
	}

	// Initialize Sample depending on how the alignment needs to be processed.
	mSWEEP::ConstructSample(alignment, args.value<size_t>("iters"), args.value<size_t>("bootstrap-count"), args.value<size_t>("seed"), bin_reads, sample);

	log.flush();
      } else {
	try {
	  // Reading both likelihoods and log_ec_counts from file.
	  log << "  reading likelihoods from file" << '\n';
	  if (reference->get_n_groupings() > 1) {
	    throw std::runtime_error("Using more than one grouping with --read-likelihood is not yet implemented.");
	  }

	  cxxio::In infile(args.value<std::string>("read-likelihood"));
	  log_likelihoods = std::make_unique<mSWEEP::LL_WOR21<double, int64_t>>();
	  log_likelihoods->from_file(reference->n_groups(i), &infile.stream());
	} catch (std::exception &e) {
	  finalize("Reading the likelihoods failed:\n  " + std::string(e.what()) + "\nexiting\n", log, true);
	  return 1;
	}
      }

      try {
	// Write the likelihood to disk here if it was requested.
	if (args.value<bool>("write-likelihood") || args.value<bool>("write-likelihood-bitseq"))
	  log_likelihoods->write((args.value<bool>("write-likelihood-bitseq") ? "bitseq" : "mSWEEP"), out.likelihoods((args.value<bool>("write-likelihood-bitseq") ? "bitseq" : "mSWEEP")));
      } catch (std::exception &e) {
	finalize("Writing the likelihood to file failed:\n  " + std::string(e.what()) + "\nexiting\n", log, true);
	return 1;
      }

      std::vector<std::string> estimated_reference_names;
      std::vector<std::string> zero_reference_names;
      // Start the abundance estimation part
      if (args.value<bool>("no-fit-model")) {
	log << "Skipping relative abundance estimation (--no-fit-model toggled)" << '\n';
      } else {
	log << "Estimating relative abundances" << '\n';

	// Prior parameters
	std::vector<double> prior_counts(log_likelihoods->log_mat().get_rows(), 1.0); // Default is all = 1.0
	if (CmdOptionPresent(argv, argv+argc, "--alphas")) {
	  if (args.value<std::vector<double>>("alphas").size() != log_likelihoods->log_mat().get_rows()) {
	    finalize("Error: --alphas must have the same number of values as there are groups.", log, true);
	    return 1;
	  }
	  prior_counts = std::move(args.value<std::vector<double>>("alphas"));
	}

	try {
	  // Run estimation
	  sample->store_probs(rcg_optl(args, log_likelihoods->log_mat(), log_likelihoods->log_counts(), prior_counts, log));
	} catch (std::exception &e) {
	  finalize("Estimating relative abundances failed:\n  " + std::string(e.what()) + "\nexiting\n", log, true);
	  return 1;
	}

	if (args.value<bool>("run-rate")) {
	    std::cerr << "WARNING: --run-rate is an experimental option that has not been thoroughly tested and is subject to change.\n" << std::endl;
	    sample->dirichlet_kld(log_likelihoods->log_counts());
	}

	if (args.value<size_t>("min-hits") > 0) {
	    std::cerr << "WARNING: --min-hits > 0 is an experimental option that has not been thoroughly tested and is subject to change.\n" << std::endl;
	}

	// Run binning if requested and write results to files.
	  // Turn the probs into relative abundances
	  if (args.value<std::string>("algorithm") == "rcgcpu") {
	      sample->store_abundances(rcgpar::mixture_components(sample->get_probs(), log_likelihoods->log_counts()));
	  } else {
	      sample->store_abundances(rcgpar::mixture_components_torch(sample->get_probs(), log_likelihoods->log_counts()));
	  }

	  if (args.value<size_t>("min-hits") > 0) {
	      for (size_t j = 0; j < reference->group_names(i).size(); ++j) {
		  if (log_likelihoods->groups_considered()[j]) {
		      estimated_reference_names.push_back(reference->group_names(i)[j]);
		  } else {
		      zero_reference_names.push_back(reference->group_names(i)[j]);
		  }
	      }
	  } else {
	      estimated_reference_names = reference->group_names(i);
	  }
	  // Bin the reads if requested
	  if (bin_reads) {
	    std::vector<std::string> target_names;
	    if (CmdOptionPresent(argv, argv+argc, "--target-groups")) {
	      target_names = std::move(args.value<std::vector<std::string>>("target-groups"));
	    } else {
	      target_names = estimated_reference_names;
	    }
	    if (CmdOptionPresent(argv, argv+argc, "--min-abundance")) {
	      mGEMS::FilterTargetGroups(estimated_reference_names, sample->get_abundances(), args.value<double>("min-abundance"), &target_names);
	    }
	    std::vector<std::vector<uint32_t>> bins;
	    try {
	      if (bootstrap_mode) {
		mSWEEP::BinningBootstrap* bs = static_cast<mSWEEP::BinningBootstrap*>(&(*sample));
		bins = std::move(mGEMS::BinFromMatrix(bs->get_aligned_reads(), sample->get_abundances(), sample->get_probs(), estimated_reference_names, &target_names));
	      } else {
		mSWEEP::BinningSample* bs = static_cast<mSWEEP::BinningSample*>(&(*sample));
		bins = std::move(mGEMS::BinFromMatrix(bs->get_aligned_reads(), sample->get_abundances(), sample->get_probs(), estimated_reference_names, &target_names));
	      }
	    } catch (std::exception &e) {
	      finalize("Binning the reads failed:\n  " + std::string(e.what()) + "\nexiting\n", log, true);
	      return 1;
	    }

	    for (size_t j = 0; j < bins.size(); ++j) {
	      try {
		mGEMS::WriteBin(bins[j], *out.bin(target_names[j]));
	      } catch (std::exception &e) {
		finalize("Writing the bin for target group " + target_names[j] + " failed:\n  " + std::string(e.what()) + "\nexiting\n", log, true);
		return 1;
	      }
	    }
	  }

	  try {
	    // Write ec_probs
	    if (args.value<bool>("print-probs")) {
	      // Note: this ignores the printing_output variable because
	      // we might want to print the probs even when writing to
	      // pipe them somewhere.
		if (args.value<size_t>("min-hits") > 0) {
		    sample->write_probs2(estimated_reference_names, zero_reference_names, &std::cout);
		} else {
		    sample->write_probs(estimated_reference_names, &std::cout);
		}
	    }
	    if (args.value<bool>("write-probs")) {
		if (args.value<size_t>("min-hits") > 0) {
		    sample->write_probs2(estimated_reference_names, zero_reference_names, out.probs());
		} else {
		    sample->write_probs(estimated_reference_names, out.probs());
		}
	    }
	  } catch (std::exception &e) {
	    finalize("Writing the probabilities failed:\n  " + std::string(e.what()) + "\nexiting\n", log, true);
	    return 1;
	  }

	// Bootstrap the ec_counts and estimate from the bootstrapped data if required
	if (bootstrap_mode) {
	  log << "Running estimation with " << args.value<size_t>("iters") << " bootstrap iterations" << '\n';
	  for (uint16_t k = 0; k < args.value<size_t>("iters"); ++k) {
	    // Bootstrap the counts
	    log << "Bootstrap" << " iter " << k + 1 << "/" << args.value<size_t>("iters") << '\n';
	    std::vector<double> resampled_counts;
	      resampled_counts = std::move(static_cast<mSWEEP::BootstrapSample*>(&(*sample))->resample_counts());

	    // Estimate with the bootstrapped counts
	    // Reuse ec_probs since it has already been processed
	    try {
	      sample->store_probs(rcg_optl(args, log_likelihoods->log_mat(), resampled_counts, prior_counts, log));
	    } catch (std::exception &e) {
	      finalize("Bootstrap iteration " + std::to_string(k) + "/" + std::to_string(args.value<size_t>("iters")) + " failed:\n  " + std::string(e.what()) + "\nexiting\n", log, true);
	      return 1;
	    }
		if (args.value<std::string>("algorithm") == "rcgcpu") {
		    sample->store_abundances(rcgpar::mixture_components(sample->get_probs(), resampled_counts));
		} else {
		    sample->store_abundances(rcgpar::mixture_components_torch(sample->get_probs(), resampled_counts));
		}
	  }
	}
      }

      // Write relative abundances
      if (!args.value<bool>("no-fit-model")) {
	try {
	  if (sample->get_rate_run()) {
	      const std::vector<double> &log_kld = sample->get_log_klds();
	      const std::vector<double> &RATE = sample->get_rates();
	      const std::vector<double> &relative_abundances = sample->get_abundances();

	      std::ostream *of = out.abundances();
	      // Write relative abundances to &of,
	      if (of->good()) {
		  (*of) << "#mSWEEP_version:" << '\t' << MSWEEP_BUILD_VERSION << '\n';
		  (*of) << "#num_reads:" << '\t' << sample->get_n_reads() << '\n';
		  (*of) << "#num_aligned:" << '\t' << sample->get_counts_total() << '\n';
		  (*of) << "#c_id" << '\t' << "mean_theta" << '\t' << "RATE" << '\t' << "KLD" << '\n';
		  size_t n_targets = estimated_reference_names.size() + zero_reference_names.size();
		  for (size_t i = 0; i < n_targets; ++i) {
		      if (i < estimated_reference_names.size()) {
			  double KLD = std::exp(log_kld[i]);
			  (*of) << estimated_reference_names[i] << '\t' << relative_abundances[i] << '\t' << RATE[i] << '\t' << KLD <<'\n';
		      } else {
			  (*of) << zero_reference_names[i - estimated_reference_names.size()] << '\t' << (double)0.0 << '\t' << (double)0.0 << '\t' << (double)0.0 <<'\n';
		      }
		  }
		  of->flush();
	      } else {
		  throw std::runtime_error("Can't write to abundances file.");
	      }
	  } else {
	    if (args.value<size_t>("min-hits") > 0) {
		sample->write_abundances2(estimated_reference_names, zero_reference_names, out.abundances());
	    } else {
		sample->write_abundances(estimated_reference_names, out.abundances());
	    }
	  }
	} catch (std::exception &e) {
	  finalize("Writing the relative abundances failed:\n  " + std::string(e.what()) + "\nexiting\n", log, true);
	  return 1;
	}
      }
      if (n_groupings > 1 && i < n_groupings - 1) {
	out.next_grouping();
      }
  }
  finalize("", log);
  return 0;
}

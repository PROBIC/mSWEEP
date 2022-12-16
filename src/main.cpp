#include <iostream>
#include <vector>
#include <memory>
#include <fstream>

#include "rcgpar.hpp"
#include "cxxargs.hpp"
#include "msweep_log.hpp"

#include "mSWEEP.hpp"
#include "Sample.hpp"
#include "Reference.hpp"

#include "Matrix.hpp"

#include "version.h"
#include "openmp_config.hpp"
#include "mpi_config.hpp"

#define MSWEEP_OPENMP_SUPPORT 1

#if defined(MSWEEP_OPENMP_SUPPORT) && (MSWEEP_OPENMP_SUPPORT) == 1
#include <omp.h>
#include <algorithm>
#pragma omp declare reduction(vec_double_plus : std::vector<double> :	\
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
#endif

void WriteResults(const cxxargs::Arguments &args, const std::unique_ptr<Sample> &sample, const Grouping &grouping, const uint16_t n_groupings, const uint16_t current_grouping) {
  cxxio::Out of;

  bool printing_output = args.value<std::string>('o').empty();
  // Set output file name correctly
  // for backwards compatibility with v1.4.0 or older
  std::string outfile = args.value<std::string>('o');
  if (n_groupings > 1 && !printing_output) {
    outfile += "_";
    outfile += std::to_string(current_grouping);
  }

  // Write likelihoods
  if (args.value<bool>("write-likelihood") || args.value<bool>("write-likelihood-bitseq")) {
    std::string ll_outfile(outfile);
    ll_outfile += (args.value<bool>("write-likelihood-bitseq") ? "_bitseq" : "");
    ll_outfile += "_likelihoods.txt";
    if (args.value<bool>("gzip-probs")) {
      ll_outfile += ".gz";
      of.open_compressed(ll_outfile);
    } else {
      of.open(ll_outfile);
    }
    if (args.value<bool>("write-likelihood-bitseq")) {
      sample->write_likelihood_bitseq(grouping.get_n_groups(), of.stream());
    } else {
      sample->write_likelihood(grouping.get_n_groups(), of.stream());
    }
  }

  // Relative abundances
  if (!args.value<bool>("no-fit-model")) {
    std::string abundances_outfile(outfile);
    if (!printing_output) {
      std::string abundances_outfile = outfile + "_abundances.txt";
      of.open(abundances_outfile);
    }
    if (args.value<size_t>("iters") > 1) {
      BootstrapSample* bs = static_cast<BootstrapSample*>(&(*sample));
      bs->write_bootstrap(grouping.get_names(), args.value<size_t>("iters"), (printing_output ? std::cout : of.stream()));
    } else {
      sample->write_abundances(grouping.get_names(), (printing_output ? std::cout : of.stream()));
    }
  }

  // Probability matrix
  if (args.value<bool>("print-probs") && !args.value<bool>("no-fit-model")) {
    sample->write_probabilities(grouping.get_names(), std::cout);
  }
  if (args.value<bool>("write-probs") && !args.value<bool>("no-fit-model")) {
    std::string probs_outfile(outfile);
    probs_outfile += "_probs.csv";
    if (args.value<bool>("gzip-probs")) {
      probs_outfile += ".gz";
      of.open_compressed(probs_outfile);
    } else {
      of.open(probs_outfile);
    }
    sample->write_probabilities(grouping.get_names(), (outfile.empty() ? std::cout : of.stream()));
  }
}

void PrintCitationInfo() {
  std::cerr << "Please cite us as:\n"
	    << "\tMäklin T, Kallonen T, David S et al. High-resolution sweep\n"
	    << "\tmetagenomics using fast probabilistic inference [version 2;\n"
	    << "\tpeer review: 2 approved]. Wellcome Open Res 2021, 5:14\n"
	    << "\t(https://doi.org/10.12688/wellcomeopenres.15639.2)" << std::endl;
}

bool CmdOptionPresent(char **begin, char **end, const std::string &option) {
  return (std::find(begin, end, option) != end);
}

void parse_args(int argc, char* argv[], cxxargs::Arguments &args) {
  // Options for outputting the probability matrix
  args.add_long_argument<bool>("write-probs", "If specified, write the read equivalence class probabilities in a .csv matrix.", false);
  args.add_long_argument<bool>("gzip-probs", "Gzip the .csv matrix output from --write-probs and the likelihoods from --write-likelihood or --write-likelihood-bitseq.", false);
  args.add_long_argument<bool>("print-probs", "Print the read equivalence class probabilities to cout.", false);

  // Write likelihood or not
  args.add_long_argument<bool>("write-likelihood", "Write the likelihood matrix to a file if -o option is specified, print to cout if -o is not.", false);
  args.add_long_argument<bool>("write-likelihood-bitseq", "Write the likelihoods in a format can be parsed by BitSeq's (https://github.com/bitseq/bitseq) functions.", false);

  // Skip running the optimizer (useful for just writing the likelihood matrix)
  args.add_long_argument<bool>("no-fit-model", "Skip fitting the model entirely. Useful if only the likelihood matrix is required.", false);

  // Alignments have been compressed with alignment-writer
  args.add_long_argument<bool>("read-compact", "Input alignments are in compact format.", false);

  // Pseudoalignment files
  args.add_long_argument<std::string>("themisto-1", "Pseudoalignment results from Themisto for the 1st strand of paired-end reads.");
  args.add_long_argument<std::string>("themisto-2", "Pseudoalignment results from Themisto for the 2nd strand of paired-end reads.");
  // Pseudoalignments are not required if reading likelihood from a file
  args.add_long_argument<std::string>("read-likelihood", "Read in a likelihood matrix that has been written to file with the --write-likelihood toggle.");
  args.set_not_required("read-likelihood");
  if (CmdOptionPresent(argv, argv+argc, "--read-likelihood") || CmdOptionPresent(argv, argv+argc, "--themisto")) {
      args.set_not_required("themisto-1");
      args.set_not_required("themisto-2");
  }
  // Separate pseudoalignment files are not required if supplied via a list
  args.add_long_argument<std::vector<std::string>>("themisto", "Single themisto alignment file.");
  args.set_not_required("themisto");

  // How to merge paired alignments
  args.add_long_argument<std::string>("themisto-mode", "How to merge Themisto pseudoalignments for paired-end reads (default: intersection).", "intersection");

  // Cluster indicators
  args.add_short_argument<std::string>('i', "Group identifiers file.");

  // Output prefix
  args.add_short_argument<std::string>('o', "Prefix for output files written from mSWEEP.", "");
  // TODO: check that the directory exists (outside of this function) using the following code:
  //   args.outfile = std::string(GetCmdOption(argv, argv+argc, "-o"));
  // if (args.outfile.find("/") != std::string::npos) {
  //   std::string outfile_dir = args.outfile;
  //   outfile_dir.erase(outfile_dir.rfind("/"), outfile_dir.size());
  //   cxxio::directory_exists(outfile_dir);
  // }

  // Number of iterations to run bootstrapping for
  args.add_long_argument<size_t>("iters", "Number of times to rerun estimation with bootstrapped alignments.", (size_t)1);

  // Seed for bootstrapping
  args.add_long_argument<size_t>("seed", "Seed for the random generator used in bootstrapping (default: random).");
  args.set_not_required("seed");

  // How many reads to resample when bootstrapping
  args.add_long_argument<size_t>("bootstrap-count", "How many reads to resample when bootstrapping (default: number of reads)");
  args.set_not_required("bootstrap-count");

  // Number of threads for parallel estimation
  args.add_short_argument<size_t>('t', "How many threads to use in abundance estimation.", (size_t)1);

  // Tolerance for abundance estimation convergence
  args.add_long_argument<double>("tol", "Optimization has converged when the bound changes less than the given tolerance.", (double)0.000001);

  // Maximum iterations to run the optimizer for
  args.add_long_argument<size_t>("max-iters", "Maximum number of iterations to run the gradient optimizer.", (size_t)5000);

  // Mean fraction of aligned sequences for the likelihood
  args.add_short_argument<double>('q', "Tuning parameter for the likelihood (mean seqs read should align against if from group).", 0.65);

  // Dispersion term for likelihood
  args.add_short_argument<double>('e', "Tuning parameter for the likelihood (dispersion term in the beta-binomial part of the likelihood).", 0.01);

  // Prior parameters for estimation
  args.add_long_argument<std::vector<double>>("alphas", "Prior parameters for the relative abundances, supply as comma-separated values (default: all 1.0).");
  args.set_not_required("alphas");

  // Print output or not.
  args.add_long_argument<bool>("verbose", "Print status messages to cerr.", false);

  args.add_long_argument<bool>("help", "Print the help message.", false);
  if (CmdOptionPresent(argv, argv+argc, "--help")) {
      std::cerr << "\n" + args.help() << '\n' << '\n';
  }
  args.parse(argc, argv);

  if (!CmdOptionPresent(argv, argv+argc, "--themisto")) {
    args.set_val<std::vector<std::string>>("themisto", std::vector<std::string>({ args.value<std::string>("themisto-1"), args.value<std::string>("themisto-2") }));
  }
}

void finalize(const std::string &msg, Log &log, bool abort = false) {
  if (abort != 0)  {
#if defined(MSWEEP_MPI_SUPPORT) && (MSWEEP_MPI_SUPPORT) == 1
    MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  }
#if defined(MSWEEP_MPI_SUPPORT) && (MSWEEP_MPI_SUPPORT) == 1
  MPI_Finalize();
#endif
  log.flush();
}

seamat::DenseMatrix<double> rcg_optl(const cxxargs::Arguments &args, const seamat::Matrix<double> &ll_mat, const std::vector<double> &log_ec_counts, const std::vector<double> &prior_counts, Log &log) {
  // Wrapper for calling rcgpar with omp or mpi depending on config
  std::ofstream of; // Silence output from ranks > 1 with an empty ofstream
#if defined(MSWEEP_MPI_SUPPORT) && (MSWEEP_MPI_SUPPORT) == 1
  // MPI parallellization
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  const seamat::DenseMatrix<double> &ec_probs = rcgpar::rcg_optl_mpi(ll_mat, log_ec_counts, prior_counts, args.value<double>("tol"), args.value<size_t>("max-iters"), (rank == 0 && args.value<bool>("verbose") ? log.stream() : of));

#else
  // OpenMP parallelllization
  const seamat::DenseMatrix<double> &ec_probs = rcgpar::rcg_optl_omp(ll_mat, log_ec_counts, prior_counts, args.value<double>("tol"), args.value<size_t>("max-iters"), (args.value<bool>("verbose") ? log.stream() : of));
#endif
  return ec_probs;
}

int main (int argc, char *argv[]) {
  int rank = 0; // If MPI is not supported then we are always on the root process
  Log log(std::cerr, CmdOptionPresent(argv, argv+argc, "--verbose"), false);

#if defined(MSWEEP_MPI_SUPPORT) && (MSWEEP_MPI_SUPPORT) == 1
  // Initialize MPI
  int rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS) {
    finalize("MPI initialization failed\n\n", log);
    return 1;
  }
  int ntasks;
  rc=MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
  rc=MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif

  // Parse command-line arguments
  cxxargs::Arguments args("mSWEEP", "Usage: mSWEEP --themisto-1 <forwardPseudoalignments> --themisto-2 <reversePseudoalignments> -i <groupIndicatorsFile>");

  log << "mSWEEP-" << MSWEEP_BUILD_VERSION << " abundance estimation" << '\n';
  if (CmdOptionPresent(argv, argv+argc, "--version")) {
    log.status(std::string("mSWEEP-") + std::string(MSWEEP_BUILD_VERSION));
  }
  if (CmdOptionPresent(argv, argv+argc, "--cite")) {
    PrintCitationInfo();
  }

  try {
    log << "Parsing arguments" << '\n';
    parse_args(argc, argv, args);
  } catch (const std::exception &e) {
    finalize("Error in parsing arguments:\n  " + std::string(e.what()) + "\nexiting\n", log);
    finalize("", log);
    return 1;
  }

  if (CmdOptionPresent(argv, argv+argc, "--version") || CmdOptionPresent(argv, argv+argc, "--help") || CmdOptionPresent(argv, argv+argc, "--cite")) {
    finalize("", log);
    return 0;
  }

#if defined(MSWEEP_OPENMP_SUPPORT) && (MSWEEP_OPENMP_SUPPORT) == 1
  omp_set_num_threads(args.value<size_t>('t'));
#endif

  Reference reference;
  log << "Reading the input files" << '\n';
  try {
    if (rank == 0) { // Only root reads in data
      log << "  reading group indicators" << '\n';
      ReadGroupIndicators(args.value<std::string>('i'), &reference);
      if (reference.get_n_groupings() > 1) {
	throw std::runtime_error("Using more than one grouping is currently unsupported.");
      }
      log << "  read " << reference.get_n_refs() << " group indicators" << '\n';
    }
  } catch (std::exception &e) {
    finalize("Reading the input files failed:\n  " + std::string(e.what()) + "\nexiting\n", log, true);
    return 1;
  }

  std::unique_ptr<Sample> sample;

  bool bootstrap_mode = args.value<size_t>("iters") > (size_t)1;
  if (bootstrap_mode) {
    sample.reset(new BootstrapSample(args.value<size_t>("seed")));
  } else {
    sample.reset(new Sample(reference));
  }
  bool likelihood_mode = CmdOptionPresent(argv, argv+argc, "--read-likelihood");

  try {
    if (!likelihood_mode) {
	log << "  reading pseudoalignments" << '\n';
	ReadPseudoalignments(args.value<std::vector<std::string>>("themisto"), args.value<std::string>("themisto-mode"), args.value<bool>("read-compact"), reference, sample);
	sample->process_aln(bootstrap_mode);
	log << "  read " << sample->num_ecs() << " unique alignments" << '\n';
	log.flush();
    } else {
      ReadLikelihoodFromFile(args.value<std::string>("read-likelihood"), reference, log.stream(), sample);
    }
  } catch (std::exception &e) {
    finalize("Reading the input files failed:\n  " + std::string(e.what()) + "\nexiting\n", log, true);
    return 1;
  }

  // Estimate abundances with all groupings that were provided
  uint16_t n_groupings;
  if (rank == 0) { // rank 0
    n_groupings = reference.get_n_groupings();
  }
#if defined(MSWEEP_MPI_SUPPORT) && (MSWEEP_MPI_SUPPORT) == 1
  // Only root process has read in the input.
  MPI_Bcast(&n_groupings, 1, MPI_UINT16_T, 0, MPI_COMM_WORLD);
#endif

  for (uint16_t i = 0; i < n_groupings; ++i) {
      // Send the number of groups from root to all processes
      uint16_t n_groups;
      if (rank == 0) // rank 0
	n_groups = reference.get_grouping(i).get_n_groups();
#if defined(MSWEEP_MPI_SUPPORT) && (MSWEEP_MPI_SUPPORT) == 1
      MPI_Bcast(&n_groups, 1, MPI_UINT16_T, 0, MPI_COMM_WORLD);
#endif

      log << "Building log-likelihood array" << '\n';
      if (rank == 0 && !likelihood_mode) // rank 0
	ConstructLikelihood(args.value<double>('q'), args.value<double>('e'), reference.get_grouping(i), reference.get_group_indicators(i), sample);

      // Process the reads accordingly
      if (args.value<bool>("no-fit-model")) {
	log << "Skipping relative abundance estimation (--no-fit-model toggled)" << '\n';
      } else {
	log << "Estimating relative abundances" << '\n';

	// Prior parameters
	std::vector<double> prior_counts(n_groups, 1.0); // Default is all = 1.0
	if (CmdOptionPresent(argv, argv+argc, "--alphas")) {
	  if (args.value<std::vector<double>>("alphas").size() != n_groups) {
	    finalize("Error: --alphas must have the same number of values as there are groups.", log, true);
	    return 1;
	  }
	  prior_counts = std::move(args.value<std::vector<double>>("alphas"));
	}

	// Run estimation
	sample->ec_probs = rcg_optl(args, sample->ll_mat, sample->log_ec_counts, prior_counts, log);
	if (rank == 0) // rank 0
	  sample->relative_abundances = rcgpar::mixture_components(sample->ec_probs, sample->log_ec_counts);

	if (bootstrap_mode) {
	  // Bootstrap the ec_counts and estimate from the bootstrapped data
	  log << "Running estimation with " << args.value<size_t>("iters") << " bootstrap iterations" << '\n';
	  BootstrapSample* bs = static_cast<BootstrapSample*>(&(*sample));
	  if (rank == 0)
	    bs->init_bootstrap();

	  for (uint16_t k = 0; k < args.value<size_t>("iters"); ++k) {
	    // Bootstrap the counts
	    std::vector<double> resampled_log_ec_counts;
	    log << "Bootstrap" << " iter " << k + 1 << "/" << args.value<size_t>("iters") << '\n';
	    if (rank == 0) {
	      size_t bootstrap_count = bs->counts_total;
	      if (CmdOptionPresent(argv, argv+argc, "--bootstrap-count"))
		bootstrap_count = args.value<size_t>("bootstrap-count");
	      resampled_log_ec_counts = bs->resample_counts(bootstrap_count);
	    }

	    // Estimate with the bootstrapped counts
	    const seamat::DenseMatrix<double> &bootstrapped_ec_probs = rcg_optl(args, bs->ll_mat, resampled_log_ec_counts, prior_counts, log);
	    if (rank == 0)
	      bs->bootstrap_results.emplace_back(rcgpar::mixture_components(bootstrapped_ec_probs, resampled_log_ec_counts));
	  }
	}
      }
      // Write results to file from the root process
      if (rank == 0)
	WriteResults(args, sample, reference.get_grouping(i), n_groupings, i);
  }
  finalize("", log);
  // sample.release();
  return 0;
}

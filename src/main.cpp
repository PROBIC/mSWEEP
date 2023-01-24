#include <iostream>
#include <vector>
#include <memory>
#include <fstream>

#include "rcgpar.hpp"
#include "cxxargs.hpp"
#include "msweep_log.hpp"
#include "bin_reads.h"

#include "mSWEEP.hpp"
#include "Sample.hpp"
#include "Reference.hpp"
#include "likelihood.hpp"

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
  // Print output or not.
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
  // TODO: check that the directory exists (outside of this function) using the following code:
  //   args.outfile = std::string(GetCmdOption(argv, argv+argc, "-o"));
  // if (args.outfile.find("/") != std::string::npos) {
  //   std::string outfile_dir = args.outfile;
  //   outfile_dir.erase(outfile_dir.rfind("/"), outfile_dir.size());
  //   cxxio::directory_exists(outfile_dir);
  // }

  // Run the mGEMS binning algorithm
  args.add_long_argument<bool>("bin-reads", "Run the mGEMS binning algorithm and write bins to the directory `-o` points to (default: false).", false);
  args.add_long_argument<std::vector<std::string>>("target-groups", "Only extract these groups, supply as comma separated list (default: extract all groups).");
  args.set_not_required("target-groups");
  args.set_not_required("target-groups");
  args.add_long_argument<double>("min-abundance", "Only extract groups that have a relative abundance higher than this value (default: 0).\n\nOutput options:");
  args.set_not_required("min-abundance");

  // Options for outputting the probability matrix
  args.add_long_argument<bool>("write-probs", "If specified, write the estimated read-to-group probabilities to a file with \"_probs.csv\" suffix (default:false).", false);
  args.add_long_argument<bool>("print-probs", "Print the read equivalence class probabilities to cout even if `-o` is given (default: false).", false);

  // Write likelihood or not
  args.add_long_argument<bool>("write-likelihood", "Write the internal likelihood matrix to a file with \"_likelihoods.txt\" suffix (default: false).", false);
  args.add_long_argument<bool>("write-likelihood-bitseq", "Write the likelihoods in a format can be parsed by BitSeq's (https://github.com/bitseq/bitseq) functions (default: false).", false);

  // Toggle compression
  args.add_long_argument<bool>("gzip-probs", "Compress the output from --write-probs, --write-likelihood, or --write-likelihood-bitseq with zlib (default: false).\n\nInput options:", false);

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
  args.add_long_argument<double>("tol", "Optimization terminates when the bound changes by less than the given tolerance (default: 0.000001).\n\nBootstrapping options:", (double)0.000001);

  // Number of iterations to run bootstrapping for
  args.add_long_argument<size_t>("iters", "Number of times to rerun estimation with bootstrapped alignments (default: 0).", (size_t)1);
  // Seed for bootstrapping
  args.add_long_argument<size_t>("seed", "Seed for the random generator used in bootstrapping (default: random).");
  args.set_not_required("seed");
  // How many reads to resample when bootstrapping
  args.add_long_argument<size_t>("bootstrap-count", "How many pseudoalignments to resample when bootstrapping (default: number of reads).\n\nLikelihood options:");
  args.set_not_required("bootstrap-count");

  // Mean fraction of aligned sequences for the likelihood
  args.add_short_argument<double>('q', "Mean for the beta-binomial component (default: 0.65).", 0.65);
  // Dispersion term for likelihood
  args.add_short_argument<double>('e', "Dispersion term for the beta-binomial component (default: 0.01).", 0.01);
  // Prior parameters for estimation
  args.add_long_argument<std::vector<double>>("alphas", "Prior counts for the relative abundances, supply as comma-separated nonzero values (default: all 1.0).");
  args.set_not_required("alphas");

  if (CmdOptionPresent(argv, argv+argc, "--help")) {
      std::cerr << "\n" + args.help() << '\n' << '\n';
  }
  args.parse(argc, argv);

  if (!CmdOptionPresent(argv, argv+argc, "--themisto") && CmdOptionPresent(argv, argv+argc, "--themisto-1") && CmdOptionPresent(argv, argv+argc, "--themisto-2")) {
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
  if (!msg.empty())
    std::cerr << msg;
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

      seamat::DenseMatrix<double> log_likelihoods;
      std::vector<double> log_ec_counts;
      try {
	if (!likelihood_mode) {
	  // TODO: handle MPI (wrap in rank == 0?)
	  log << "  reading pseudoalignments" << '\n';
	  const telescope::GroupedAlignment &alignment = ReadPseudoalignments(args.value<std::vector<std::string>>("themisto"), args.value<std::string>("themisto-mode"), reference);
	  sample->process_aln(alignment, bootstrap_mode);

	  log_ec_counts.resize(alignment.n_ecs(), 0);
#pragma omp parallel for schedule(static)
	  for (uint32_t i = 0; i < alignment.n_ecs(); ++i) {
	    log_ec_counts[i] = std::log(alignment.reads_in_ec(i));
	  }

	  log << "  read " << sample->num_ecs() << " unique alignments" << '\n';
	  log.flush();

	  log << "Building log-likelihood array" << '\n';
	  if (rank == 0 && !likelihood_mode) // rank 0 TODO: this could be built distributed to avoid communicating it to other processses.
	    log_likelihoods = std::move(likelihood_array_mat(alignment, reference.get_grouping(i), args.value<double>('q'), args.value<double>('e')));
	  log.flush();
	} else {
	  log_likelihoods = std::move(ReadLikelihoodFromFile(args.value<std::string>("read-likelihood"), reference, log.stream(), &log_ec_counts));
	}

	if (args.value<bool>("write-likelihood") || args.value<bool>("write-likelihood-bitseq") && rank == 0) {
	  cxxio::Out of;
	  std::string ll_outfile(args.value<std::string>('o'));
	  ll_outfile += (args.value<bool>("write-likelihood-bitseq") ? "_bitseq" : "");
	  ll_outfile += "_likelihoods.txt";
	  if (args.value<bool>("gzip-probs")) {
	    ll_outfile += ".gz";
	    of.open_compressed(ll_outfile);
	  } else {
	    of.open(ll_outfile);
	  }
	  if (args.value<bool>("write-likelihood-bitseq")) {
	    WriteLikelihoodBitSeq(log_likelihoods, log_ec_counts, reference.get_grouping(i).get_n_groups(), of.stream());
	  } else {
	    WriteLikelihood(log_likelihoods, log_ec_counts, reference.get_grouping(i).get_n_groups(), of.stream());
	  }
	  of.close();
	}
      } catch (std::exception &e) {
	finalize("Reading the input files failed:\n  " + std::string(e.what()) + "\nexiting\n", log, true);
	return 1;
      }

      std::vector<double> relative_abundances;
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
	const seamat::DenseMatrix<double> &ec_probs = rcg_optl(args, log_likelihoods, log_ec_counts, prior_counts, log);
	if (rank == 0) { // rank 0
	  relative_abundances = std::move(rcgpar::mixture_components(ec_probs, log_ec_counts));

	  // Bin the reads if requested
	  if (args.value<bool>("bin-reads")) {
	    std::vector<std::string> target_names;
	    if (CmdOptionPresent(argv, argv+argc, "--target-groups")) {
	      target_names = std::move(args.value<std::vector<std::string>>("target-groups"));
	    } else {
	      target_names = reference.get_grouping(i).get_names();
	    }
	    if (CmdOptionPresent(argv, argv+argc, "--min-abundance")) {
	      mGEMS::FilterTargetGroups(reference.get_grouping(i).get_names(), relative_abundances, args.value<double>("min-abundance"), &target_names);
	    }
	    const std::vector<std::vector<uint32_t>> &bins = mGEMS::BinFromMatrix(sample->pseudos, relative_abundances, ec_probs, reference.get_grouping(i).get_names(), &target_names);
	    std::string outfile_dir = args.value<std::string>('o');
	    outfile_dir.erase(outfile_dir.rfind("/"), outfile_dir.size()); // TODO check that path contains a /
	    for (size_t j = 0; j < bins.size(); ++j) {
	      cxxio::Out of(outfile_dir + '/' + target_names[j] + ".bin");
	      mGEMS::WriteBin(bins[j], of.stream());
	    }
	  }

	  bool printing_output = args.value<std::string>('o').empty();
	  // Write the results
	  std::string outfile = args.value<std::string>('o');
	  if (n_groupings > 1 && !printing_output) {
	    outfile += "_";
	    outfile += std::to_string(i);
	  }

	  cxxio::Out of;
	  // Write relative abundances
	  if (!args.value<bool>("no-fit-model")) {
	    std::string abundances_outfile(outfile);
	    if (!printing_output) {
	      std::string abundances_outfile = outfile + "_abundances.txt";
	      of.open(abundances_outfile);
	    }
	    if (args.value<size_t>("iters") > 1) {
	      // TODO this needs to be written at the end after bootstrapping is done. (store relative abundances in BootstrapSample)
	      //BootstrapSample* bs = static_cast<BootstrapSample*>(&(*sample));
	      //bs->write_bootstrap(reference.get_grouping(i).get_names(), args.value<size_t>("iters"), (printing_output ? std::cout : of.stream()));
	    } else {
	      WriteAbundances(relative_abundances, reference.get_grouping(i).get_names(), sample->counts_total, (printing_output ? std::cout : of.stream()));
	    }
	  }

	  // Write ec_probs
	  if (args.value<bool>("print-probs") && !args.value<bool>("no-fit-model")) {
	    WriteProbabilities(ec_probs, reference.get_grouping(i).get_names(), std::cout);
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
	    WriteProbabilities(ec_probs, reference.get_grouping(i).get_names(), (outfile.empty() ? std::cout : of.stream()));
	  }
	  of.close();
	}

	if (bootstrap_mode) {
	  // Bootstrap the ec_counts and estimate from the bootstrapped data
	  log << "Running estimation with " << args.value<size_t>("iters") << " bootstrap iterations" << '\n';
	  BootstrapSample* bs = static_cast<BootstrapSample*>(&(*sample));
	  if (rank == 0)
	    bs->init_bootstrap();

	  for (uint16_t k = 0; k < args.value<size_t>("iters"); ++k) {
	    // Bootstrap the counts
	    log << "Bootstrap" << " iter " << k + 1 << "/" << args.value<size_t>("iters") << '\n';
	    if (rank == 0) {
	      size_t bootstrap_count = bs->counts_total;
	      if (CmdOptionPresent(argv, argv+argc, "--bootstrap-count"))
		bootstrap_count = args.value<size_t>("bootstrap-count");
	      log_ec_counts = std::move(bs->resample_counts(bootstrap_count));
	    }

	    // Estimate with the bootstrapped counts
	    const seamat::DenseMatrix<double> &bootstrapped_ec_probs = rcg_optl(args, log_likelihoods, log_ec_counts, prior_counts, log);
	    if (rank == 0)
	      bs->bootstrap_results.emplace_back(rcgpar::mixture_components(bootstrapped_ec_probs, log_ec_counts));
	  }
	  cxxio::Out of;
	  std::string outfile = args.value<std::string>('o');
	  std::string abundances_outfile(outfile);
	  if (!args.value<std::string>('o').empty()) {
	    std::string abundances_outfile = outfile + "_abundances.txt";
	    of.open(abundances_outfile);
	  }
	  WriteBootstrappedAbundances(relative_abundances, bs->bootstrap_results, reference.get_grouping(i).get_names(), args.value<size_t>("iters"), bs->counts_total, (args.value<std::string>('o').empty() ? std::cout : of.stream()));
	  of.close();
	}
      }
  }
  finalize("", log);
  // sample.release();
  return 0;
}

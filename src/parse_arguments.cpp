#include "parse_arguments.hpp"

#include <algorithm>
#include <iostream>
#include <exception>

#include "cxxio.hpp"

void PrintHelpMessage() {
  std::cerr << "Usage: mSWEEP -f <pseudomappingFile> -i <clusterIndicators> [OPTIONS]\n"
	    << "Estimates the group abundances in a sample.\n\n"
	    << "Options:\n"
	    << "\t--themisto-1 <themistoPseudoalignment1>\n"
	    << "\tPseudoalignment results from Themisto for the 1st strand of paired-end reads.\n"
	    << "\t--themisto-2 <themistoPseudoalignment2>\n"
	    << "\tPseudoalignment results from Themisto for the 2nd strand of paired-end reads.\n"
	    << "\n"
	    << "\t-i <clusterIndicators>\n"
	    << "\tGroup identifiers file. Must be supplied.\n"
	    << "\t-o <outputFile>\n"
	    << "\tOutput file (folder when estimating from a batch) to write results in.\n"
    	    << "\t-t <nrThreads>\n"
	    << "\tHow many threads to use. (default: 1)\n"
	    << "\n"
	    << "\t--themisto-mode <PairedEndMergeMode>\n"
	    << "\tHow to merge Themisto pseudoalignments for paired-end reads	(default: intersection).\n"
    	    << "\t--themisto-index <ThemistoIndex>\n"
	    << "\tPath to the Themisto index the pseudoalignment was performed against (optional).\n"
	    << "\n"
    	    << "\t--fasta <ReferenceSequences>\n"
	    << "\tPath to the reference sequences the pseudoalignment index was constructed from (optional, Themisto v1.2.0 or older only)\n"
    	    << "\t--groups-list <groupIndicatorsList>\n"
	    << "\tTable containing names of the reference sequences (1st column) and their group assignments (2nd column) (optional)\n"
    	    << "\t--groups-delimiter <groupIndicatorsListDelimiter>\n"
	    << "\tDelimiter character for the --groups option (optional, default: tab)\n"
	    << "\n"
	    << "\t--iters <nrIterations>\n"
	    << "\tNumber of times to rerun estimation with bootstrapped alignments (default: 1)\n"
    	    << "\t--bootstrap-count <nrBootstrapCount>\n"
	    << "\tHow many reads to resample when bootstrapping (integer, default: all)\n"
    	    << "\t--seed <BootstrapSeed>\n"
	    << "\tSeed for the random generator used in bootstrapping (default: random)\n"    
	    << "\n"
            << "\t--write-probs\n"
            << "\tIf specified, write the read equivalence class probabilities in a .csv matrix\n"
            << "\t--print-probs\n"
            << "\tPrint the read equivalence class probabilities to cout\n"
            << "\t--write-likelihood\n"
            << "\tWrite the likelihood matrix to a file if -o option is specified, print to cout if -o is not.\n"
            << "\t--write-likelihood-bitseq\n"
            << "\tWrite the likelihoods in a format can be parsed by BitSeq's (https://github.com/bitseq/bitseq) functions.\n"
            << "\t--gzip-probs\n"
            << "\tGzip the .csv matrix output from --write-probs and the likelihoods from --write-likelihood or --write-likelihood-bitseq.\n"
	    << "\n"
	    << "\t--read-likelihood\n"
	    << "\tRead in a likelihood matrix that has been written to file with the --write-likelihood toggle.\n"
	    << "\n"
            << "\t--read-compact\n"
            << "\tInput alignments are in compact format.\n"
	    << "\t--help\n"
	    << "\tPrint this message.\n"
	    << "\t--version\n"
	    << "\tPrint the version number.\n"
	    << "\t--cite\n"
	    << "\tPrint citation information.\n"
	    << "\n\tELBO optimization and modeling\n"
	    << "\t--no-fit-model\n"
	    << "\tSkip fitting the model entirely. Useful if only the likelihood matrix is required.\n"
	    << "\t--tol <tolerance>\n"
	    << "\tOptimization has converged when the bound changes less than the given tolerance.\n"
	    << "\t--max-iters\n"
	    << "\tMaximum number of iterations to run the gradient optimizer.\n"
	    << "\t-q <meanFraction>\n"
	    << "\tFraction of the sequences in a group that the mean is set to."
	    << " (default: 0.65)\n"
	    << "\t-e <dispersionTerm>\n"
	    << "\tCalibration term in the likelihood function."
	    << " (default: 0.01)" << std::endl;
}

void PrintCitationInfo() {
  std::cerr << "Please cite us as:\n"
	    << "\tMäklin T, Kallonen T, David S et al. High-resolution sweep\n"
	    << "\tmetagenomics using fast probabilistic inference [version 2;\n"
	    << "\tpeer review: 2 approved]. Wellcome Open Res 2021, 5:14\n"
	    << "\t(https://doi.org/10.12688/wellcomeopenres.15639.2)" << std::endl;
}

char* GetCmdOption(char **begin, char **end, const std::string &option) {
  char **it = std::find(begin, end, option);
  return ((it != end && ++it != end) ? *it : 0);
}

bool CmdOptionPresent(char **begin, char **end, const std::string &option) {
  return (std::find(begin, end, option) != end);
}

double ParseDoubleOption(char **begin, char **end, const std::string &option) {
  double opt;
  char* char_opt = GetCmdOption(begin, end, option);
  if (char_opt != 0) {
    opt = std::stod(std::string(char_opt));
  } else {
    throw std::runtime_error(option + " specified but no value given");
  }
  return(opt);
}

void ParseArguments(int argc, char *argv[], Arguments &args) {
  if (argc < 3) {
    throw std::runtime_error("Error: Specify at least the infiles and the indicators file.\n");
  }

  args.optimizer.write_probs = CmdOptionPresent(argv, argv+argc, "--write-probs");
  args.optimizer.gzip_probs = CmdOptionPresent(argv, argv+argc, "--gzip-probs");
  args.optimizer.print_probs = CmdOptionPresent(argv, argv+argc, "--print-probs");

  args.optimizer.write_likelihood = CmdOptionPresent(argv, argv+argc, "--write-likelihood");
  args.optimizer.write_likelihood_bitseq = CmdOptionPresent(argv, argv+argc, "--write-likelihood-bitseq");
  args.optimizer.no_fit_model = CmdOptionPresent(argv, argv+argc, "--no-fit-model");

  args.compact_alignments = CmdOptionPresent(argv, argv+argc, "--read-compact");

  if (CmdOptionPresent(argv, argv+argc, "--themisto-1") && CmdOptionPresent(argv, argv+argc, "--themisto-2")) {
    args.tinfile1 = std::string(GetCmdOption(argv, argv+argc, "--themisto-1"));
    args.tinfile2 = std::string(GetCmdOption(argv, argv+argc, "--themisto-2"));
    if (CmdOptionPresent(argv, argv+argc, "--themisto-mode")) {
      args.themisto_merge_mode = std::string(GetCmdOption(argv, argv+argc, "--themisto-mode"));
    } else {
      args.themisto_merge_mode = std::string("intersection");
    }
    if (CmdOptionPresent(argv, argv+argc, "--themisto-index")) {
      args.themisto_index_path = std::string(GetCmdOption(argv, argv+argc, "--themisto-index"));
      cxxio::directory_exists(args.themisto_index_path);
    }
  } else if (CmdOptionPresent(argv, argv+argc, "--read-likelihood")) {
      args.likelihood_file = std::string(GetCmdOption(argv, argv+argc, "--read-likelihood"));
      args.read_likelihood_mode = true;
      if (CmdOptionPresent(argv, argv+argc, "--themisto-index")) {
	args.themisto_index_path = std::string(GetCmdOption(argv, argv+argc, "--themisto-index"));
	cxxio::directory_exists(args.themisto_index_path);
      }
    } else {
    throw std::runtime_error("infile not found.");
  }

  if (CmdOptionPresent(argv, argv+argc, "-i")) {
    args.indicators_file = std::string(GetCmdOption(argv, argv+argc, "-i"));
  } else if (CmdOptionPresent(argv, argv+argc, "--indicators")) {
    args.indicators_file = std::string(GetCmdOption(argv, argv+argc, "--indicators"));
  } else if (!CmdOptionPresent(argv, argv+argc, "--fasta") || !CmdOptionPresent(argv, argv+argc, "--groups-list")) {
    throw std::runtime_error("group indicator file not found.");
  }

  if (CmdOptionPresent(argv, argv+argc, "-o")) {
    args.outfile = std::string(GetCmdOption(argv, argv+argc, "-o"));
    if (args.outfile.find("/") != std::string::npos) {
      std::string outfile_dir = args.outfile;
      outfile_dir.erase(outfile_dir.rfind("/"), outfile_dir.size());
      cxxio::directory_exists(outfile_dir);
    }
  }

  if (CmdOptionPresent(argv, argv+argc, "--iters")) {
    signed nr_iters_given = std::stoi(std::string(GetCmdOption(argv, argv+argc, "--iters")));
    args.bootstrap_mode = true;
    if (nr_iters_given < 1) {
      throw std::runtime_error("number of iterations must be greater or equal to 1");
    } else {
      args.iters = nr_iters_given;
    }
  }

  if (CmdOptionPresent(argv, argv+argc, "--seed")) {
    int32_t seed_given = std::stoi(std::string(GetCmdOption(argv, argv+argc, "--seed")));
    if (seed_given < 1) {
      throw std::runtime_error("seed must be greater or equal to 1");
    } else {
      args.seed = seed_given;
    }
  }

  if (CmdOptionPresent(argv, argv+argc, "--fasta") || CmdOptionPresent(argv, argv+argc, "--groups-list")) {
    if ((!CmdOptionPresent(argv, argv+argc, "--fasta") || !CmdOptionPresent(argv, argv+argc, "--groups-list"))) {
      throw std::runtime_error("--fasta and --groups-list must both be specified if either is present.");
    }
    args.fasta_file = std::string(GetCmdOption(argv, argv+argc, "--fasta"));
    args.groups_list_file = std::string(GetCmdOption(argv, argv+argc, "--groups-list"));
  }

  if (CmdOptionPresent(argv, argv+argc, "--groups-delimiter")) {
    std::string groups_list_delimiter = std::string(GetCmdOption(argv, argv+argc, "--groups-delimiter"));
    if (groups_list_delimiter.size() > 1) {
      throw std::runtime_error("--groups-delimiter must be a single character");
    } else {
      args.groups_list_delimiter = groups_list_delimiter.at(0);
    }
  }

  if (CmdOptionPresent(argv, argv+argc, "--bootstrap-count")) {
    signed bootstrap_count_given = std::stoi(std::string(GetCmdOption(argv, argv+argc, "--bootstrap-count")));
    if (bootstrap_count_given < 1) {
      throw std::runtime_error("bootstrap-count must be greater or equal to 1");
    } else {
      args.bootstrap_count = bootstrap_count_given;
    }
  }

  if (CmdOptionPresent(argv, argv+argc, "-t")) {
    signed nr_threads_given = std::stoi(std::string(GetCmdOption(argv, argv+argc, "-t")));
    if (nr_threads_given < 1) {
      throw std::runtime_error("number of threads must be strictly positive");
    } else {
      args.optimizer.nr_threads = nr_threads_given;
    }
  } else {
    args.optimizer.nr_threads = 1;
  }

  if (CmdOptionPresent(argv, argv+argc, "--tol")) {
    double tolerance = ParseDoubleOption(argv, argv+argc, "--tol");
    if (tolerance <= 0) {
      throw std::runtime_error("tolerance can't be negative or zero");
    } else {
      args.optimizer.tolerance = tolerance;
    }
  }

  if (CmdOptionPresent(argv, argv+argc, "--max-iters")) {
    int32_t max_iters = std::stoi(std::string(GetCmdOption(argv, argv+argc, "--max-iters")));
    if (max_iters < 0 || max_iters > 65536) {
      throw std::runtime_error("--max-iters must be at least 1 and less than 65537.");
    } else {
      args.optimizer.max_iters = (uint16_t)max_iters;
    }
  }

  if (CmdOptionPresent(argv, argv+argc, "-q")) {
    double frac_mu = ParseDoubleOption(argv, argv+argc, "-q");
    if (frac_mu <= 0.5 || frac_mu >= 1.0) {
      throw std::runtime_error("-q must be between 0.5 and 1.");
    } else {
      args.optimizer.bb_constants[0] = frac_mu;
    }
  }
  
  if (CmdOptionPresent(argv, argv+argc, "-e")) {
    double epsilon = ParseDoubleOption(argv, argv+argc, "-e");
    if (epsilon <= 0.0 || epsilon >= 2.0*(args.optimizer.bb_constants[0]) - 1.0) {
      throw std::runtime_error("-e must be greater than 0, and less than 2*q - 1");
    } else {
      args.optimizer.bb_constants[1] = epsilon;
    }
  }
}

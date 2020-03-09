#include <algorithm>
#include <iostream>
#include <exception>

#include "parse_arguments.hpp"

void PrintHelpMessage() {
  std::cerr << "Usage: mSWEEP -f <pseudomappingFile> -i <clusterIndicators> [OPTIONS]\n"
	    << "Estimates the group abundances in a sample.\n\n"
	    << "Options:\n"
	    << "\t--themisto-1 <themistoPseudoalignment1>\n"
	    << "\tPseudoalignment results from Themisto for the 1st strand of paired-end reads.\n"
	    << "\t--themisto-2 <themistoPseudoalignment2>\n"
	    << "\tPseudoalignment results from Themisto for the 2nd strand of paired-end reads.\n"
	    << "\n"
	    << "\t-f <pseudomappingFile>\n"
	    << "\tPseudoalignment output file location from kallisto. Can't be used when -b is specified.\n"
    	    << "\t-b <pseudomappingBatch>\n"
	    << "\tThe kallisto batch matrix file location. Can't be used when -f is specified.\n"
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
            << "\tPrint the equivalence class probabilities rather than writing when using --write-probs\n"    
            << "\t--gzip-probs\n"
            << "\tGzip the .csv matrix output from --write-probs\n"
	    << "\t--help\n"
	    << "\tPrint this message.\n"
	    << "\n\tELBO optimization and modeling (these seldom need to be changed)\n"
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
  if (CmdOptionPresent(argv, argv+argc, "--help")) {
    throw std::invalid_argument("");
  }
  if (argc < 3) {
    throw std::runtime_error("Error: Specify at least the infile and indicators file.\n");
  }
  std::cerr << "Parsing arguments" << std::endl;

  args.optimizer.write_probs = CmdOptionPresent(argv, argv+argc, "--write-probs");
  args.optimizer.gzip_probs = CmdOptionPresent(argv, argv+argc, "--gzip-probs");
  args.optimizer.print_probs = CmdOptionPresent(argv, argv+argc, "--print-probs");
  
  if ((CmdOptionPresent(argv, argv+argc, "-f") || CmdOptionPresent(argv, argv+argc, "--file"))  && CmdOptionPresent(argv, argv+argc, "-b")) {
    throw std::runtime_error("infile and batchfile found, specify only one");
  } else if (CmdOptionPresent(argv, argv+argc, "-f")) {
    args.infile = std::string(GetCmdOption(argv, argv+argc, "-f"));
  } else if (CmdOptionPresent(argv, argv+argc, "--file")) {
    args.infile = std::string(GetCmdOption(argv, argv+argc, "--file"));
  } else if (CmdOptionPresent(argv, argv+argc, "-b")) {
    args.batch_infile = std::string(GetCmdOption(argv, argv+argc, "-b"));
    args.batch_mode = true;
  } else if (CmdOptionPresent(argv, argv+argc, "--themisto-1") && CmdOptionPresent(argv, argv+argc, "--themisto-2")) {
    args.tinfile1 = std::string(GetCmdOption(argv, argv+argc, "--themisto-1"));
    args.tinfile2 = std::string(GetCmdOption(argv, argv+argc, "--themisto-2"));
    args.themisto_mode = true;
    if (CmdOptionPresent(argv, argv+argc, "--themisto-mode")) {
      args.themisto_merge_mode = std::string(GetCmdOption(argv, argv+argc, "--themisto-mode"));
    } else {
      args.themisto_merge_mode = std::string("intersection");
    }
  } else {
    throw std::runtime_error("infile not found.");
  }

  // Fill the kallisto_files vector
  args.kallisto_files = std::vector<std::string>((args.batch_mode ? 4 : 3));
  if (args.batch_mode) {
    args.infiles = KallistoFiles(args.batch_infile, args.batch_mode);
    args.kallisto_files[0] = args.batch_infile + "/run_info.json";
    args.kallisto_files[1] = args.batch_infile + "/matrix.ec";
    args.kallisto_files[2] = args.batch_infile + "/matrix.tsv";
    args.kallisto_files[3] = args.batch_infile + "/matrix.cells";
  } else if (!args.themisto_mode) {
    args.infiles = KallistoFiles(args.infile, args.batch_mode);
    args.kallisto_files[0] = args.infile + "/run_info.json";
    args.kallisto_files[1] = args.infile + "/pseudoalignments.ec";
    args.kallisto_files[2] = args.infile + "/pseudoalignments.tsv";
  }
  if (CmdOptionPresent(argv, argv+argc, "--compressed-input") && !args.themisto_mode) {
    for (size_t i = 1; i < args.kallisto_files.size(); ++i) {
      args.kallisto_files[i] += ".gz";
    }
  }

  if (CmdOptionPresent(argv, argv+argc, "-i")) {
    args.indicators_file = std::string(GetCmdOption(argv, argv+argc, "-i"));
  } else if (CmdOptionPresent(argv, argv+argc, "--indicators")) {
    args.indicators_file = std::string(GetCmdOption(argv, argv+argc, "--indicators"));
  } else {
    throw std::runtime_error("group indicator file not found.");
  }

  if (CmdOptionPresent(argv, argv+argc, "-o")) {
    args.outfile = std::string(GetCmdOption(argv, argv+argc, "-o"));
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
    unsigned max_iters = std::stoi(std::string(GetCmdOption(argv, argv+argc, "--max-iters")));
    if (max_iters < 0) {
      throw std::runtime_error("--max-iters must be at least 1");
    } else {
      args.optimizer.max_iters = max_iters;
    }
  }

  if (CmdOptionPresent(argv, argv+argc, "-q")) {
    double frac_mu = ParseDoubleOption(argv, argv+argc, "-q");
    if (frac_mu <= 0.5 || frac_mu >= 1.0) {
      throw std::runtime_error("-q must be between 0.5 and 1.");
    } else {
      args.params[0] = frac_mu;
    }
  }
  
  if (CmdOptionPresent(argv, argv+argc, "-e")) {
    double epsilon = ParseDoubleOption(argv, argv+argc, "-e");
    if (epsilon <= 0.0 || epsilon >= 2.0*(args.params[0]) - 1.0) {
      throw std::runtime_error("-e must be greater than 0, and less than 2*q - 1");
    } else {
      args.params[1] = epsilon;
    }
  }
}

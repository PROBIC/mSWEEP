#include <algorithm>
#include <iostream>
#include <exception>

#include "parse_arguments.hpp"

void PrintHelpMessage() {
  std::cerr << "Usage: mSWEEP -f <pseudomappingFile> -i <clusterIndicators> [OPTIONS]\n"
	    << "Estimates the group abundances in a sample.\n\n"
	    << "Options:\n"
	    << "\t-f <pseudomappingFile>\n"
	    << "\tThe pseudoalignment output file location. Can't be used when -b is specified.\n"
    	    << "\t-b <pseudomappingBatch>\n"
	    << "\tThe kallisto batch matrix file location. Can't be used when -f is specified.\n"
	    << "\t-i <clusterIndicators>\n"
	    << "\tGroup identifiers file. Must be supplied.\n"
	    << "\t-o <outputFile>\n"
	    << "\tOutput file (folder when estimating from a batch) to write results in.\n"
    	    << "\t-t <nrThreads>\n"
	    << "\tHow many threads to use when processing a batch matrix (default: 1)\n"
            << "\n\t--write-probs\n"
            << "\tIf specified, write the read equivalence class probabilities in a .csv matrix\n"
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
  
  if ((CmdOptionPresent(argv, argv+argc, "-f") || CmdOptionPresent(argv, argv+argc, "--file"))  && CmdOptionPresent(argv, argv+argc, "-b")) {
    throw std::runtime_error("infile and batchfile found, specify only one");
  } else if (CmdOptionPresent(argv, argv+argc, "-f")) {
    args.infile = std::string(GetCmdOption(argv, argv+argc, "-f"));
  } else if (CmdOptionPresent(argv, argv+argc, "--file")) {
    args.infile = std::string(GetCmdOption(argv, argv+argc, "--file"));
  } else if (CmdOptionPresent(argv, argv+argc, "-b")) {
    args.batch_infile = std::string(GetCmdOption(argv, argv+argc, "-b"));
  } else {
    throw std::runtime_error("infile not found.");
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

  if (CmdOptionPresent(argv, argv+argc, "-t")) {
    if (!CmdOptionPresent(argv, argv+argc, "-b")) {
      std::cerr << "  parallel processing a single sample is not supported" << std::endl;
      args.nr_threads = 1;
    } else {
      signed nr_threads_given = std::stoi(std::string(GetCmdOption(argv, argv+argc, "-t")));
      if (nr_threads_given < 1) {
	throw std::runtime_error("number of threads must be strictly positive");
      } else {
	args.nr_threads = nr_threads_given;
      }
    }
  } else {
    args.nr_threads = 1;
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

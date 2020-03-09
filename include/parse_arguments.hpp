#ifndef MSWEEP_PARSE_ARGUMENTS_HPP
#define MSWEEP_PARSE_ARGUMENTS_HPP

#include <vector>
#include <string>

#include "KallistoFiles.hpp"

struct OptimizerArgs {
  uint16_t max_iters = 5000;
  double tolerance = 1e-06;

  std::vector<double> alphas;
  bool write_probs;
  bool gzip_probs;
  bool print_probs;
  unsigned nr_threads = 1;
};

struct Arguments {
  OptimizerArgs optimizer;

  KallistoFiles infiles;
  std::string infile;
  std::string batch_infile;
  std::string indicators_file;
  std::string outfile;
  std::string tinfile1;
  std::string tinfile2;
  std::vector<std::string> kallisto_files;

  bool batch_mode = false;
  bool bootstrap_mode = false;
  bool compressed_input = false;
  bool themisto_mode = false;

  std::string themisto_merge_mode = "union";

  unsigned iters = 1;
  uint32_t bootstrap_count = 0;
  int32_t seed = -1;
  double params[2] = { 0.65, 0.01 };

  // 0 = single sample, 1 = batch input, 2 = bootstrap single sample
  uint8_t run_mode() { return 0 | (this->bootstrap_mode << 1) | (this->batch_mode << 0); };
};

void ParseArguments(int argc, char *argv[], Arguments &args);
void PrintHelpMessage();

#endif

#ifndef PARSE_ARGUMENTS_H
#define PARSE_ARGUMENTS_H

#include <vector>
#include <string>

struct OptimizerArgs {
  unsigned max_iters = 5000;
  double tolerance = 1e-06;

  std::vector<double> alphas;
  bool write_probs;
};

struct Arguments {
  OptimizerArgs optimizer;

  std::string infile;
  std::string batch_infile;
  std::string indicators_file;
  std::string outfile;
  std::vector<std::string> kallisto_files;

  bool batch_mode = false;
  bool bootstrap_mode = false;

  unsigned nr_threads = 1;
  unsigned iters = 1;
  double params[2] = { 0.65, 0.01 };

  // 0 = single sample, 1 = batch input, 2 = bootstrap single sample
  uint8_t run_mode() { return 0 | (this->bootstrap_mode << 1) | (this->batch_mode << 0); };
};

void ParseArguments(int argc, char *argv[], Arguments &args);
void PrintHelpMessage();

#endif

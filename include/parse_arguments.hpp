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

  unsigned nr_threads = 1;
  double params[2] = { 0.65, 0.01 };
};

void ParseArguments(int argc, char *argv[], Arguments &args);
void PrintHelpMessage();

#endif

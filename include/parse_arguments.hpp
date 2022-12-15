#ifndef MSWEEP_PARSE_ARGUMENTS_HPP
#define MSWEEP_PARSE_ARGUMENTS_HPP

#include <vector>
#include <string>

struct OptimizerArgs {
  uint16_t max_iters = 5000;
  double tolerance = 1e-06;
  double bb_constants[2] = { 0.65, 0.01 };

  std::vector<double> alphas;
  bool write_probs;
  bool write_likelihood;
  bool write_likelihood_bitseq;
  bool no_fit_model;
  bool gzip_probs;
  bool print_probs;
  unsigned nr_threads = 1;
};

struct Arguments {
  OptimizerArgs optimizer;

  std::string indicators_file;
  std::string outfile;
  std::string tinfile1;
  std::string tinfile2;
  std::string likelihood_file;

  std::string fasta_file;
  std::string groups_list_file;
  char groups_list_delimiter = '\t';

  bool bootstrap_mode = false;
  bool read_likelihood_mode = false;
  bool compact_alignments = false;

  std::string themisto_merge_mode = "union";

  uint16_t iters = 1;
  uint32_t bootstrap_count = 0;
  int32_t seed = -1;
};

void ParseArguments(int argc, char *argv[], Arguments &args);
void PrintHelpMessage();
void PrintCitationInfo();
bool CmdOptionPresent(char **begin, char **end, const std::string &option);

#endif

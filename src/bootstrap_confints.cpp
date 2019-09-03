#include <iostream>
#include <fstream>

#include "bootstrap_confints.hpp"
#include "process_reads.hpp"
#include "version.h"

void write_bootstrap(const std::vector<std::string> &cluster_indicators_to_string, const std::vector<std::vector<double>> &abundances, std::string &outfile, unsigned iters, unsigned counts_total) {
  // Write relative abundances to a file,
  // outputs to std::cout if outfile is empty.
  std::streambuf *buf;
  std::ofstream of;
  if (outfile.empty()) {
    buf = std::cout.rdbuf();
  } else {
    outfile += "_abundances.txt";
    of.open(outfile);
    buf = of.rdbuf();
  }
  std::ostream out(buf);
  out << "#mSWEEP_version:" << '\t' << _BUILD_VERSION << '\n';
  out << "#total_hits:" << '\t' << counts_total << '\n';
  out << "#bootstrap_iters:" << '\t' << iters << '\n';
  out << "#c_id" << '\t' << "mean_theta" << '\t' << "bootstrap_mean_thetas" << '\n';

  for (size_t i = 0; i < cluster_indicators_to_string.size(); ++i) {
    out << cluster_indicators_to_string[i] << '\t';
    for (unsigned j = 0; j <= iters; ++j) {
      out << abundances.at(j).at(i) << (j == iters ? '\n' : '\t');
    }
  }
  out << std::endl;
  if (!outfile.empty()) {
    of.close();
  }
}

void bootstrap_abundances(const Arguments &args, std::string &outfile, std::mt19937_64 &gen, Reference &reference, Sample &sample, ThreadPool &pool) {
  std::vector<std::future<std::vector<double>>> abus;
  // Init the bootstrap variables
  std::cerr << "Building log-likelihood array" << std::endl;
  sample.init_bootstrap(reference.grouping);
  std::cout << "Running estimation on the original sample" << std::endl;
  for (unsigned i = 0; i <= args.iters; ++i) {
    if (i != 0) {
      std::cout << "Running bootstrap iteration " << i << " out of " << args.iters << std::endl;
    }
    // Run the estimation multiple times without writing anything
    abus.emplace_back(pool.enqueue(&ProcessBootstrap, reference, sample, sample.ec_counts, args.optimizer));
    // Resample the pseudoalignment counts (here because we want to include the original)
    sample.resample_counts(gen);
  }
  std::vector<std::vector<double>> results(args.iters + 1);
  for (size_t i = 0; i <= args.iters; ++i) {
    results[i] = abus.at(i).get();
  }

  pool.enqueue(&write_bootstrap, reference.group_names, results, outfile, args.iters, sample.total_counts());
}

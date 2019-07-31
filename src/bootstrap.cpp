#include <fstream>
#include <iostream>
#include <future>
#include <random>

#include "bootstrap.hpp"
#include "Reference.hpp"
#include "process_reads.hpp"
#include "parse_arguments.hpp"
#include "thread_pool.hpp"

void write_bootstrap(const std::vector<std::string> &cluster_indicators_to_string, const std::vector<std::vector<double>> &abundances, std::string &outfile, unsigned iters) {
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
  out << "#c_id" << '\t' << "abundances" << '\t' << "bootstrap_abundances" << '\n';

  for (size_t i = 0; i < cluster_indicators_to_string.size(); ++i) {
    out << cluster_indicators_to_string[i] << '\t';
    for (unsigned j = 0; j < iters; ++j) {
      out << abundances.at(j).at(i) << (j == iters - 1 ? '\n' : '\t');
    }
  }
  out << std::endl;
  if (!outfile.empty()) {
    of.close();
  }
}

std::unordered_map<std::string, std::vector<std::vector<double>>> bootstrap_abundances(const std::vector<Sample> &bitfields, Reference &reference, ThreadPool &pool, Arguments &args) {
    std::unordered_map<std::string, std::vector<std::vector<double>>> results;
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::cerr << "Running estimation with " << args.iters << " bootstrap iterations" << '\n';
    for (auto bitfield : bitfields) {
      // Store results in this
      std::vector<std::future<std::vector<double>>> abus;
      // Init the bootstrap variables
      std::cerr << "Building log-likelihood array" << std::endl;
      bitfield.init_bootstrap(reference.grouping);
      for (unsigned i = 0; i < args.iters; ++i) {
	// Run the estimation multiple times without writing anything
	abus.emplace_back(pool.enqueue(&BootstrapIter, reference, bitfield, bitfield.ec_counts, args.optimizer));
	// Resample the pseudoalignment counts (here because we want to include the original)
	bitfield.resample_counts(gen);
      }
      std::string name = (args.batch_mode ? bitfield.cell_name() : "0");
      results.insert(std::make_pair(name, std::vector<std::vector<double>>()));
      for (unsigned i = 0; i < args.iters; ++i) {
	results.at(name).emplace_back(abus.at(i).get());
      }
    }
    return results;
}

#include "openmp_config.hpp"

#include <fstream>
#include <iostream>
#include <future>
#include <random>

#include "Reference.hpp"
#include "parse_arguments.hpp"
#include "thread_pool.hpp"
#include "rcg.hpp"
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

std::vector<double> bootstrap_iter(Reference &reference, Sample &sample, std::vector<double> ec_counts, OptimizerArgs args) {
  // Process pseudoalignments but return the abundances rather than writing.
  std::cerr << "Estimating relative abundances" << std::endl;
  sample.ec_probs = rcg_optl_mat(sample.ll_mat, sample.total_counts(), ec_counts, args.alphas, args.tolerance, args.max_iters);
  const std::vector<double> &abundances = sample.group_abundances();

  return abundances;
}

//std::unordered_map<std::string, std::vector<std::vector<double>>> bootstrap_abundances(const std::vector<Sample> &bitfields, Reference &reference, ThreadPool &pool, Arguments &args) {
BootstrapResults bootstrap_abundances(const std::vector<Sample> &bitfields, Reference &reference, ThreadPool &pool, Arguments &args) {
  //    std::unordered_map<std::string, std::vector<std::vector<double>>> results;
  BootstrapResults results;
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::cerr << "Running estimation with " << args.iters << " bootstrap iterations" << '\n';
    for (auto bitfield : bitfields) {
      // Which sample are we processing?
      std::string name = (args.batch_mode ? bitfield.cell_name() : "0");
      std::cout << "Processing " << (args.batch_mode ? name : "the sample") << std::endl;
      // Store results in this
#if defined(MSWEEP_OPENMP_SUPPORT) && (MSWEEP_OPENMP_SUPPORT) == 0
      std::vector<std::future<std::vector<double>>> abus;
#else
      std::vector<std::vector<double>> abus;
#endif
      // Init the bootstrap variables
      std::cerr << "Building log-likelihood array" << std::endl;
      bitfield.init_bootstrap(reference.grouping);
      for (unsigned i = 0; i <= args.iters; ++i) {
	if (i > 0) {
	  std::cout << "Bootstrapping alignment counts" << '\n';
	  std::cout << "  iter: " << i << "/" << args.iters << std::endl;
	}
	// Run the estimation multiple times without writing anything
#if defined(MSWEEP_OPENMP_SUPPORT) && (MSWEEP_OPENMP_SUPPORT) == 0
	abus.emplace_back(pool.enqueue(&bootstrap_iter, reference, bitfield, bitfield.ec_counts, args.optimizer));
#else
	abus.emplace_back(bootstrap_iter(reference, bitfield, bitfield.ec_counts, args.optimizer));
#endif
	// Resample the pseudoalignment counts (here because we want to include the original)
	bitfield.resample_counts(gen);
      }
      //      results.insert(std::make_pair(name, std::vector<std::vector<double>>()));
      results.insert(name, bitfield.total_counts(), std::vector<std::vector<double>>());
      for (unsigned i = 0; i <= args.iters; ++i) {
	//	results.at(name).emplace_back(abus.at(i).get());
#if defined(MSWEEP_OPENMP_SUPPORT) && (MSWEEP_OPENMP_SUPPORT) == 0
	results.insert_iter(name, abus.at(i).get());
#else
	results.insert_iter(name, abus.at(i));
#endif
      }
    }
    return results;
}

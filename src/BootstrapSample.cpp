#include "Sample.hpp"

#include "likelihood.hpp"
#include "rcg.hpp"
#include "version.h"

#include "rcgpar.hpp"
#include "bxzstr.hpp"
#include "cxxio.hpp"

BootstrapSample::BootstrapSample(const int32_t seed) {
  if (seed == -1) {
    std::random_device rd;
    this->gen = std::mt19937_64(rd());
  } else {
    this->gen = std::mt19937_64(seed);
  }
}

std::vector<double> BootstrapSample::resample_counts(const uint32_t how_many, std::mt19937_64 &generator) {
  std::vector<uint32_t> tmp_counts(num_ecs());
  for (uint32_t i = 0; i < how_many; ++i) {
    uint32_t ec_id = ec_distribution(generator);
    tmp_counts[ec_id] += 1;
  }
  std::vector<double> resampled_log_ec_counts(num_ecs());
#pragma omp parallel for schedule(static)
  for (uint32_t i = 0; i < num_ecs(); ++i) {
    resampled_log_ec_counts[i] = std::log(tmp_counts[i]);
  }
  return resampled_log_ec_counts;
}

void BootstrapSample::bootstrap_iter(const std::vector<double> &resampled_log_ec_counts,
				     const std::vector<double> &alpha0, const double tolerance,
				     const uint16_t max_iters) {
  // Process pseudoalignments but return the abundances rather than writing.
  const rcgpar::Matrix<double> &bootstrapped_ec_probs = rcgpar::rcg_optl_omp(this->ll_mat, resampled_log_ec_counts, alpha0, tolerance, max_iters, std::cerr);
  this->relative_abundances.emplace_back(rcgpar::mixture_components(bootstrapped_ec_probs, resampled_log_ec_counts));
}

void BootstrapSample::bootstrap_abundances(const Grouping &grouping, const Arguments &args) {
  std::cerr << "Running estimation with " << args.iters << " bootstrap iterations" << '\n';
  // Which sample are we processing?
  std::string name = (args.batch_mode ? cell_name() : "0");
  std::cout << "Processing " << (args.batch_mode ? name : "the sample") << std::endl;

  // Initialize ec_distribution for bootstrapping
  ec_distribution = std::discrete_distribution<uint32_t>(pseudos.ec_counts.begin(), pseudos.ec_counts.end());

  // Clear the abundances in case we're estimating the same sample again.
  this->relative_abundances = std::vector<std::vector<double>>();

  // Estimate the relative abundances from the input data
  std::cerr << "Estimating relative abundances without bootstrapping" << std::endl;
  this->ec_probs = rcgpar::rcg_optl_omp(this->ll_mat, this->log_ec_counts, args.optimizer.alphas, args.optimizer.tolerance, args.optimizer.max_iters, std::cerr);
  this->relative_abundances.emplace_back(rcgpar::mixture_components(this->ec_probs, this->log_ec_counts));

  for (uint16_t i = 0; i <= args.iters; ++i) {
    std::cout << "Bootstrap" << " iter " << i << "/" << args.iters << std::endl;

    // Resample the pseudoalignment counts
    const std::vector<double> resampled_log_ec_counts = resample_counts((args.bootstrap_count == 0 ? this->get_counts_total() : args.bootstrap_count), gen);

    // Estimate with the resampled counts
    bootstrap_iter(resampled_log_ec_counts, args.optimizer.alphas, args.optimizer.tolerance, args.optimizer.max_iters);
  }
}


void BootstrapSample::write_bootstrap(const std::vector<std::string> &cluster_indicators_to_string, std::string outfile, const unsigned iters, const bool batch_mode) const {
  // Write relative abundances to a file,
  // outputs to std::cout if outfile is empty.
  outfile = (outfile.empty() || !batch_mode ? outfile : outfile + '/' + cell_name());
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
  out << "#mSWEEP_version:" << '\t' << MSWEEP_BUILD_VERSION << '\n';
  out << "#total_hits:" << '\t' << this->get_counts_total() << '\n';
  out << "#bootstrap_iters:" << '\t' << iters << '\n';
  out << "#c_id" << '\t' << "mean_theta" << '\t' << "bootstrap_mean_thetas" << '\n';

  for (size_t i = 0; i < cluster_indicators_to_string.size(); ++i) {
    out << cluster_indicators_to_string[i] << '\t';
    for (unsigned j = 0; j <= iters; ++j) {
      out << relative_abundances[j][i] << (j == iters ? '\n' : '\t');
    }
  }
  if (!outfile.empty()) {
    of.close();
  }
}

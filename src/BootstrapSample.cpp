#include "Sample.hpp"

#include "rcgpar.hpp"

#include "version.h"

BootstrapSample::BootstrapSample(const int32_t seed) {
  if (seed == -1) {
    std::random_device rd;
    this->gen = std::mt19937_64(rd());
  } else {
    this->gen = std::mt19937_64(seed);
  }
}

std::vector<double> BootstrapSample::resample_counts(const uint32_t how_many) {
  std::vector<uint32_t> tmp_counts(num_ecs());
  for (uint32_t i = 0; i < how_many; ++i) {
    uint32_t ec_id = ec_distribution(this->gen);
    tmp_counts[ec_id] += 1;
  }
  std::vector<double> resampled_log_ec_counts(num_ecs());
#pragma omp parallel for schedule(static)
  for (uint32_t i = 0; i < num_ecs(); ++i) {
    resampled_log_ec_counts[i] = std::log(tmp_counts[i]);
  }
  return resampled_log_ec_counts;
}

void BootstrapSample::init_bootstrap() {
  // Clear the bootstrap abundances in case we're estimating the same sample again.
  this->bootstrap_results = std::vector<std::vector<double>>();

  // Initialize ec_distribution for bootstrapping
  ec_distribution = std::discrete_distribution<uint32_t>(pseudos.ec_counts.begin(), pseudos.ec_counts.end());
}

void BootstrapSample::write_bootstrap(const std::vector<std::string> &cluster_indicators_to_string,
				      const uint16_t iters, std::ostream &of) const {
  // Write relative abundances to a file,
  // outputs to std::cout if outfile is empty.
  if (of.good()) {
    of << "#mSWEEP_version:" << '\t' << MSWEEP_BUILD_VERSION << '\n';
    of << "#total_hits:" << '\t' << this->get_counts_total() << '\n';
    of << "#bootstrap_iters:" << '\t' << iters << '\n';
    of << "#c_id" << '\t' << "mean_theta" << '\t' << "bootstrap_mean_thetas" << '\n';

    for (size_t i = 0; i < cluster_indicators_to_string.size(); ++i) {
      of << cluster_indicators_to_string[i] << '\t';
      of << relative_abundances[i] << '\t';
      for (uint16_t j = 0; j < iters; ++j) {
	of << bootstrap_results[j][i] << (j == iters - 1 ? '\n' : '\t');
      }
    }
    of.flush();
  } else {
    throw std::runtime_error("Could not write to abundances file.");
  }
}

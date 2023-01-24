#include "Sample.hpp"

#include "rcgpar.hpp"

#include "version.h"

void BootstrapSample::init_bootstrap(const telescope::GroupedAlignment &alignment) {
  // Clear the bootstrap abundances in case we're estimating the same sample again.
  this->bootstrap_results = std::vector<std::vector<double>>();

  // Initialize ec_distribution for bootstrapping
  std::vector<uint32_t> weights(this->num_ecs);
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < this->num_ecs; ++i) {
    weights[i] = alignment.reads_in_ec(i);
  }
  ec_distribution = std::discrete_distribution<uint32_t>(weights.begin(), weights.end());
}

BootstrapSample::BootstrapSample(const telescope::GroupedAlignment &alignment, const int32_t seed) {
  if (seed == -1) {
    std::random_device rd;
    this->gen = std::mt19937_64(rd());
  } else {
    this->gen = std::mt19937_64(seed);
  }
  this->num_ecs = alignment.n_ecs();
  this->init_bootstrap(alignment);
}

std::vector<double> BootstrapSample::resample_counts(const uint32_t how_many) {
  std::vector<uint32_t> tmp_counts(this->num_ecs);
  for (uint32_t i = 0; i < how_many; ++i) {
    uint32_t ec_id = ec_distribution(this->gen);
    tmp_counts[ec_id] += 1;
  }
  std::vector<double> resampled_log_ec_counts(this->num_ecs);
#pragma omp parallel for schedule(static)
  for (uint32_t i = 0; i < this->num_ecs; ++i) {
    resampled_log_ec_counts[i] = std::log(tmp_counts[i]);
  }
  return resampled_log_ec_counts;
}

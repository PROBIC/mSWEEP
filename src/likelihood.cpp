#include "likelihood.hpp"

#include <vector>
#include <cmath>

inline double lbeta(double x, double y) {
  return(std::lgamma(x) + std::lgamma(y) - std::lgamma(x + y));
}

inline double log_bin_coeff(uint16_t n, uint16_t k) {
  return (std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1));
}

inline double ldbb_scaled(uint16_t k, uint16_t n, double alpha, double beta) {
  return (log_bin_coeff(n, k) + lbeta(k + alpha, n - k + beta) - lbeta(n + alpha, beta));
}

void precalc_lls(const Grouping &grouping, Sample *sample) {
  uint16_t max_size = 0;
  for (uint32_t i = 0; i < grouping.n_groups; ++i) {
    max_size = (grouping.sizes[i] > max_size ? grouping.sizes[i] : max_size);
  }

  //  Matrix<double> lls(grouping.n_groups, max_size + 1, -4.60517); // log(0.01) = -4.60517
  sample->ll_mat.mat.resize(grouping.n_groups, std::vector<double>(max_size + 1, -4.60517));
  sample->ll_mat.rows = grouping.n_groups;
  sample->ll_mat.cols = max_size + 1;
#pragma omp parallel for schedule(static)
  for (uint32_t i = 0; i < grouping.n_groups; ++i) {
    for (uint16_t j = 1; j <= max_size; ++j) {
      sample->ll_mat.mat[i][j] = ldbb_scaled(j, grouping.sizes[i], grouping.bb_params[i][0], grouping.bb_params[i][1]) - 0.01005034; // log(0.99) = -0.01005034
    }
  }
}

void likelihood_array_mat(Sample &sample, const Grouping &grouping) {
  uint32_t num_ecs = sample.num_ecs();
  precalc_lls(grouping, &sample);

  sample.counts.resize(grouping.n_groups, std::vector<uint16_t>(sample.num_ecs(), 0));
#pragma omp parallel for schedule(static)
  for (uint32_t j = 0; j < num_ecs; ++j) {
    const std::vector<uint16_t> &counts = sample.group_counts(grouping.indicators, j, grouping.n_groups);
    for (uint32_t i = 0; i < grouping.n_groups; ++i) {
      sample.counts[i][j] = counts[i];
    }
  }
  sample.clear_configs();
}

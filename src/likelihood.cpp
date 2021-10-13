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

void precalc_lls(const Grouping &grouping, const double bb_constants[2], Matrix<double> *ll_mat) {
  const std::vector<std::array<double, 2>> &bb_params = grouping.bb_parameters(bb_constants);

  uint16_t max_size = 0;
  for (uint32_t i = 0; i < grouping.n_groups; ++i) {
    max_size = (grouping.sizes[i] > max_size ? grouping.sizes[i] : max_size);
  }

  ll_mat->resize(grouping.n_groups, max_size + 1, -4.60517);
#pragma omp parallel for schedule(static)
  for (uint32_t i = 0; i < grouping.n_groups; ++i) {
    for (uint16_t j = 1; j <= max_size; ++j) {
      (*ll_mat)(i, j) = ldbb_scaled(j, grouping.sizes[i], bb_params[i][0], bb_params[i][1]) - 0.01005034; // log(0.99) = -0.01005034
    }
  }
}

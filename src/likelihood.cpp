#include "likelihood.hpp"

#include "rcgpar.hpp"
#include "Matrix.hpp"

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

void precalc_lls(const Grouping &grouping, const double bb_constants[2], seamat::DenseMatrix<double> &ll_mat) {
  const std::vector<std::array<double, 2>> &bb_params = grouping.bb_parameters(bb_constants);
  uint32_t n_groups = grouping.get_n_groups();

  uint16_t max_size = 0;
  for (uint32_t i = 0; i < n_groups; ++i) {
    max_size = (grouping.get_sizes()[i] > max_size ? grouping.get_sizes()[i] : max_size);
  }

  ll_mat.resize(n_groups, max_size + 1, -4.60517);
#pragma omp parallel for schedule(static) shared(ll_mat)
  for (uint32_t i = 0; i < n_groups; ++i) {
    for (uint16_t j = 1; j <= max_size; ++j) {
      ll_mat(i, j) = ldbb_scaled(j, grouping.get_sizes()[i], bb_params[i][0], bb_params[i][1]) - 0.01005034; // log(0.99) = -0.01005034
    }
  }
}

seamat::DenseMatrix<double> likelihood_array_mat(const telescope::GroupedAlignment &pseudos, const Grouping &grouping, const double tol, const double frac_mu) {
  uint32_t num_ecs = pseudos.n_ecs();
  uint16_t n_groups = grouping.get_n_groups();

  seamat::DenseMatrix<double> precalc_lls_mat;
  double bb_constants[2] = { tol, frac_mu };
  precalc_lls(grouping, bb_constants, precalc_lls_mat);

  seamat::DenseMatrix<double> log_likelihoods(n_groups, num_ecs, -4.60517);
#pragma omp parallel for schedule(static) shared(precalc_lls_mat)
  for (size_t j = 0; j < num_ecs; ++j) {
    for (size_t i = 0; i < n_groups; ++i) {
      log_likelihoods(i, j) = precalc_lls_mat(i, pseudos.get_group_count(i, j));
    }
  }
  return log_likelihoods;
}

#include "likelihood.hpp"

#include <vector>
#include <cmath>

inline double lbeta(double x, double y) {
  return(std::lgamma(x) + std::lgamma(y) - std::lgamma(x + y));
}

inline double log_bin_coeff(unsigned short n, unsigned short k) {
  return (std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1));
}

inline double ldbb_scaled(unsigned short k, unsigned short n, double alpha, double beta) {
  return (log_bin_coeff(n, k) + lbeta(k + alpha, n - k + beta) - lbeta(n + alpha, beta));
}

Matrix<double> precalc_lls(const Grouping &grouping) {
  unsigned short max_size = 0;
  for (unsigned short i = 0; i < grouping.n_groups; ++i) {
    max_size = (grouping.sizes[i] > max_size ? grouping.sizes[i] : max_size);
  }

  Matrix<double> lls(grouping.n_groups, max_size + 1, -4.60517); // log(0.01) = -4.60517
#pragma omp parallel for schedule(static)
  for (unsigned short i = 0; i < grouping.n_groups; ++i) {
    for (unsigned short j = 1; j <= max_size; ++j) {
      lls(i, j) = ldbb_scaled(j, grouping.sizes[i], grouping.bb_params[i][0], grouping.bb_params[i][1]) - 0.01005034; // log(0.99) = -0.01005034
    }
  }

  return lls;
}

Matrix<double> likelihood_array_mat(Sample &sample, const Grouping &grouping) {
  unsigned num_ecs = sample.num_ecs();
  const Matrix<double> &my_lls = precalc_lls(grouping);
  Matrix<double> log_likelihoods(grouping.n_groups, num_ecs, -4.60517);

  //#pragma omp parallel for schedule(static)
  // for (unsigned j = 0; j < num_ecs; ++j) {
  //   const std::vector<short unsigned> &counts = sample.group_counts(grouping.indicators, j);
  //   for (unsigned short i = 0; i < grouping.n_groups; ++i) {
  //     log_likelihoods(i, j) = my_lls(i, counts[i]);
  //   }
  // }
  sample.counts.resize(grouping.n_groups, std::vector<short unsigned>(sample.m_num_ecs, 0));
#pragma omp parallel for schedule(static)
  for (unsigned j = 0; j < num_ecs; ++j) {
    const std::vector<short unsigned> &counts = sample.group_counts(grouping.indicators, j);
    for (unsigned i = 0; i < grouping.n_groups; ++i) {
      sample.counts[i][j] = counts[i];
    }
  }

  return my_lls;
  //  return(log_likelihoods);
}

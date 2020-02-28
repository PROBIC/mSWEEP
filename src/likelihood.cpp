#include "likelihood.hpp"

#include <vector>
#include <cmath>

inline double lbeta(double x, double y) {
  return(std::lgamma(x) + std::lgamma(y) - std::lgamma(x + y));
}

inline double log_bin_coeff(unsigned n, unsigned k) {
  return (std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1));
}

inline double ldbb_scaled(unsigned k, unsigned n, double alpha, double beta) {
  return (log_bin_coeff(n, k) + lbeta(k + alpha, n - k + beta) - lbeta(n + alpha, beta));
}

Matrix<double> likelihood_array_mat(const Sample &sample, Grouping grouping) {
  // Pass grouping by copy since multiple threads accessing the values
  // (which they often do) by reference slow things down somewhat.
  unsigned num_ecs = sample.num_ecs();
  Matrix<double> log_likelihoods(grouping.n_groups, num_ecs, 0.0);

#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < num_ecs; ++i) {
    std::vector<unsigned> read_hitcounts = sample.group_counts(grouping.indicators, grouping.n_groups, i);
    for (size_t j = 0; j < grouping.n_groups; ++j) {
      unsigned cluster_hits = read_hitcounts[j];
      if (cluster_hits == 0) {
	log_likelihoods(j, i) = -4.60517; //log(0.01)
      } else if (grouping.sizes[j] == 1) {
	log_likelihoods(j, i) = 0.0;
      } else {
	log_likelihoods(j, i) = ldbb_scaled(cluster_hits, grouping.sizes[j], grouping.bb_params[j][0], grouping.bb_params[j][1]);
      }
    }
  }
  return(log_likelihoods);
}

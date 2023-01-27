// mSWEEP: Estimate abundances of reference lineages in DNA sequencing reads.
//
// MIT License
//
// Copyright (c) 2023 Probabilistic Inference and Computational Biology group @ UH
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
#ifndef MSWEEP_LIKELIHOOD_HPP
#define MSWEEP_LIKELIHOOD_HPP

#include "Matrix.hpp"
#include "telescope.hpp"

#include "Grouping.hpp"

#include <cmath>
#include <cstddef>
#include <vector>
#include <array>

#include "openmp_config.hpp"

#include "Grouping.hpp"

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

#endif

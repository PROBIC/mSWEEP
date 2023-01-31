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

#include "mSWEEP_openmp_config.hpp"

#include <cmath>
#include <cstddef>
#include <vector>
#include <array>
#include <sstream>

#include "openmp_config.hpp"

#include "Grouping.hpp"

template <typename T>
T lbeta(T x, T y) {
  return(std::lgamma(x) + std::lgamma(y) - std::lgamma(x + y));
}

template <typename T, typename V>
T log_bin_coeff(V n, V k) {
  return (std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1));
}

template <typename T, typename V>
T ldbb_scaled(V k, V n, T alpha, T beta) {
  return (log_bin_coeff<T, V>(n, k) + lbeta<T>(k + alpha, n - k + beta) - lbeta<T>(n + alpha, beta));
}

template <typename T>
class Likelihood {
public:
  virtual T& operator()(const size_t row, const size_t col) =0;
  virtual const T& operator()(const size_t row, const size_t col) const =0;

  virtual void from_file(const size_t n_groups, std::istream* infile) =0;

  virtual void write(const std::string &format, std::ostream *of) const =0;

  // Get the matrix
  virtual const seamat::Matrix<T>& log_mat() const =0;

  // Get the ec counts
  virtual const std::vector<T>& log_counts() const =0;

};

template <typename T, typename V>
class LL_WOR21 : public Likelihood<T> {
private:
  seamat::DenseMatrix<T> log_likelihoods;
  std::vector<T> log_ec_counts;
  std::vector<std::array<T, 2>> bb_params;

  seamat::DenseMatrix<T> precalc_lls(const std::vector<V> &group_sizes, const size_t n_groups) {
    V max_size = 0; // Storing the grouping can take a lot less space if it can be done with uint16_t or uint8_t.
    for (size_t i = 0; i < n_groups; ++i) {
      max_size = (group_sizes[i] > max_size ? group_sizes[i] : max_size);
    }

    seamat::DenseMatrix<T> ll_mat(n_groups, max_size + 1, -4.60517);
#pragma omp parallel for schedule(static) shared(ll_mat)
    for (size_t i = 0; i < n_groups; ++i) {
      for (V j = 1; j <= max_size; ++j) {
	ll_mat(i, j) = ldbb_scaled(j, group_sizes[i], this->bb_params[i][0], this->bb_params[i][1]) - 0.01005034; // log(0.99) = -0.01005034
      }
    }

    return ll_mat;
  }

  void fill_ll_mat(const telescope::Alignment &alignment, const std::vector<V> &group_sizes, const size_t n_groups) {
    size_t num_ecs = alignment.n_ecs();

    const seamat::DenseMatrix<T> &precalc_lls_mat = this->precalc_lls(group_sizes, n_groups);

    this->log_likelihoods.resize(n_groups, num_ecs, -4.60517); // -4.60517 = log(0.01)
#pragma omp parallel for schedule(static) shared(precalc_lls_mat)
    for (size_t j = 0; j < num_ecs; ++j) {
      for (size_t i = 0; i < n_groups; ++i) {
	this->log_likelihoods(i, j) = precalc_lls_mat(i, alignment(i, j));
      }
    }
  }

  void fill_ec_counts(const telescope::Alignment &alignment) {
    // Fill log ec counts.
    this->log_ec_counts.resize(alignment.n_ecs(), 0);
#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < alignment.n_ecs(); ++i) {
      this->log_ec_counts[i] = std::log(alignment.reads_in_ec(i));
    }
  }

  // Calculates the group-specific likelihood parameters that depend on the group sizes
  void update_bb_parameters(const std::vector<V> &group_sizes, const size_t n_groups, const double bb_constants[2]) {
    this->bb_params = std::vector<std::array<double, 2>>(n_groups);
    for (size_t i = 0; i < n_groups; ++i) {
      double e = group_sizes[i]*bb_constants[0];
      double phi = 1.0/(group_sizes[i] - e + bb_constants[1]);
      double beta = phi*(group_sizes[i] - e);
      double alpha = (e*beta)/(group_sizes[i] - e);
      this->bb_params[i] = std::array<double, 2>{ { alpha, beta } };
    }
  }

public:
  LL_WOR21() = default;

  LL_WOR21(const std::vector<V> &group_sizes, const size_t n_groups, const T tol, const T frac_mu) {
    T bb_constants[2] = { tol, frac_mu };
    this->update_bb_parameters(group_sizes, n_groups, bb_constants);
  }

  void from_grouped_alignment(const telescope::Alignment &alignment, const std::vector<V> &group_sizes, const size_t n_groups) {
    this->fill_ll_mat(alignment, group_sizes, n_groups);
    this->fill_ec_counts(alignment);
  }

  void from_file(const size_t n_groups, std::istream* infile) {
    // Have to read the likelihoods into a temporary because num_ecs is not known
    std::vector<std::vector<double>> likelihoods(n_groups, std::vector<double>());
    this->log_ec_counts.resize(0);

    if (infile->good()) {
      std::string newline;
      uint32_t line_nr = 0;
      while (std::getline(*infile, newline)) {
	++line_nr;
	std::string part;
	std::stringstream partition(newline);
	bool ec_count_col = true;
	uint32_t group_id = 0;
	while (std::getline(partition, part, '\t')) {
	  if (ec_count_col) {
	    uint32_t ec_count = std::stol(part);
	    this->log_ec_counts.emplace_back(std::log(ec_count));
	    ec_count_col = false;
	  } else {
	    likelihoods[group_id].emplace_back(std::stod(part));
	    ++group_id;
	  }
	}
      }
    } else {
      throw std::runtime_error("Could not read from the likelihoods file.");
    }
    this->log_likelihoods = std::move(likelihoods);
  }

  void write_likelihood_mSWEEP(std::ostream *of) const {
    // Write likelihoods to a file
    size_t n_ecs = this->log_ec_counts.size();
    size_t n_groups = this->log_likelihoods.get_rows();

    if (of->good()) {
      for (uint32_t i = 0; i < n_ecs; ++i){
	uint32_t ec_hit_count = std::round(std::exp(log_ec_counts[i]));
	*of << ec_hit_count << '\t';
	for (uint32_t j = 0; j < n_groups; ++j) {
	  *of << this->log_likelihoods(j, i);
	  *of << (j == n_groups - 1 ? '\n' : '\t');
	}
      }
      of->flush();
    } else {
      throw std::runtime_error("Can't write to likelihoods file.");
    }
  }

  void write_likelihood_BitSeq(std::ostream *of) const {
    // Write likelihoods to a file
    // *Note*: will write in BitSeq format!
    // Use Sample::write_likelihoods if tab-separated matrix format is needed.

    auto sum_of_exps = [](size_t accumulator, const double &val) {
      return accumulator + std::exp(val);
    };

    size_t n_ecs = this->log_ec_counts.size();
    size_t n_groups = this->log_likelihoods.get_rows();
    if (of->good()) {
      size_t counts_total = std::accumulate(this->log_ec_counts.begin(), this->log_ec_counts.end(), 0, sum_of_exps);
      *of << "# Ntotal " << counts_total << '\n';
      *of << "# Nmap " << counts_total << '\n';
      *of << "# M " << n_groups << '\n';
      *of << "# LOGFORMAT (probabilities saved on log scale.)" << '\n';
      *of << "# r_name num_alignments (tr_id prob )^*{num_alignments}" << '\n';

      uint32_t read_id = 1;
      for (uint32_t i = 0; i < n_ecs; ++i) {
	uint32_t ec_hit_count = std::round(std::exp(log_ec_counts[i]));
	for (uint32_t k = 0; k < ec_hit_count; ++k) {
	  *of << read_id << ' ';
	  *of << n_groups + 1 << ' ';
	  for (uint32_t j = 0; j < n_groups; ++j) {
	    *of << j + 1 << ' ' << this->log_likelihoods(j, i) << ' ';
	  }
	  *of << 0 << ' ' << "-10000.00" << '\n';
	  ++read_id;
	}
      }
      of->flush();
    } else {
      throw std::runtime_error("Can't write to likelihoods file (bitseq format).");
    }
  }

  void write(const std::string &format, std::ostream *of) const override {
    if (format == "bitseq") {
      this->write_likelihood_BitSeq(of);
    } else {
      this->write_likelihood_mSWEEP(of);
    }
  }

  T& operator()(const size_t row, const size_t col) override { return this->log_likelihoods(row, col); }
  const T& operator()(const size_t row, const size_t col) const override { return this->log_likelihoods(row, col); }

  // Get the matrix
  const seamat::Matrix<T>& log_mat() const override { return this->log_likelihoods; };

  // Get the ec counts
  const std::vector<T>& log_counts() const override { return this->log_ec_counts; };

};

#endif

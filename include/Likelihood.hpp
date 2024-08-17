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

#include "mSWEEP_openmp_config.hpp"

#include <cmath>
#include <cstddef>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <sstream>
#include <exception>
#include <algorithm>
#include <memory>

#include "mSWEEP_alignment.hpp"
#include "Grouping.hpp"

namespace mSWEEP {
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

  // Get vector indicating which groups were included
  virtual const std::vector<bool>& groups_considered() const =0;
};

template <typename T, typename V>
class LL_WOR21 : public Likelihood<T> {
private:
  seamat::DenseMatrix<T> log_likelihoods;
  std::vector<T> log_ec_counts;
  std::vector<std::array<T, 2>> bb_params;
  std::vector<bool> groups_mask;
  T zero_inflation;
  T bb_constants[2];

  seamat::DenseMatrix<T> precalc_lls(const std::vector<V> &group_sizes, const size_t n_groups) {
    V max_size = 0; // Storing the grouping can take a lot less space if it can be done with uint16_t or uint8_t.
    for (size_t i = 0; i < n_groups; ++i) {
      max_size = (group_sizes[i] > max_size ? group_sizes[i] : max_size);
    }

    seamat::DenseMatrix<T> ll_mat(n_groups, max_size + 1, std::log(this->zero_inflation));
#pragma omp parallel for schedule(static) shared(ll_mat)
    for (size_t i = 0; i < n_groups; ++i) {
      for (V j = 1; j <= max_size; ++j) {
	ll_mat(i, j) = ldbb_scaled(j, group_sizes[i], this->bb_params[i][0], this->bb_params[i][1]) + std::log1p(-this->zero_inflation);
      }
    }

    return ll_mat;
  }

  void fill_ll_mat(const mSWEEP::Alignment &alignment, const std::vector<V> &group_sizes, const size_t n_groups, const size_t min_hits) {
    size_t num_ecs = alignment.n_ecs();
    size_t n_targets = alignment.get_n_targets();

    size_t n_threads = 1;
#if defined(MSWEEP_OPENMP_SUPPORT) && (MSWEEP_OPENMP_SUPPORT) == 1
#pragma omp parallel
    {
	n_threads = omp_get_num_threads();
    }
#endif

    // This double loop is currently the slowest part in the input reading
    std::vector<size_t> group_counts(num_ecs*n_groups, 0);
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_ecs; ++i) {
	for (size_t j = 0; j < n_targets; ++j) {
	    if (alignment(i, j)) {
		++group_counts[alignment.get_groups()[j]*num_ecs + i];
	    }
	}
    }

    bool mask_groups = min_hits > 0;
    this->groups_mask = std::vector<bool>(n_groups, !mask_groups);
    std::vector<V> masked_group_sizes;
    std::vector<size_t> groups_pos(n_groups, 0);
    size_t n_masked_groups = 0;
    if (mask_groups) {
	std::vector<size_t> group_hit_counts(n_groups, (size_t)0);
	// Create mask identifying groups that have at least 1 alignment
#pragma omp parallel for schedule(static) reduction(vec_size_t_plus:group_hit_counts)
	for (size_t j = 0; j < n_groups; ++j) {
	    for (size_t i = 0; i < num_ecs; ++i) {
		group_hit_counts[j] += (group_counts[j*num_ecs + i] > 0) * alignment.reads_in_ec(i);
	    }
	}

	for (size_t i = 0; i < n_groups; ++i) {
	    this->groups_mask[i] = groups_mask[i] || (group_hit_counts[i] >= min_hits);
	    if (this->groups_mask[i]) {
		groups_pos[i] = n_masked_groups;
		masked_group_sizes.push_back(group_sizes[i]);
		++n_masked_groups;
	    }
	}
    } else {
	masked_group_sizes = group_sizes;
#pragma omp parallel for schedule(static)
	for (size_t i = 0; i < n_groups; ++i) {
	    groups_pos[i] = i;
	}
    }
    n_masked_groups = masked_group_sizes.size();

    this->update_bb_parameters(masked_group_sizes, n_masked_groups, this->bb_constants);
    const seamat::DenseMatrix<T> &precalc_lls_mat = this->precalc_lls(masked_group_sizes, n_masked_groups);

    this->log_likelihoods.resize(n_masked_groups, num_ecs, std::log(this->zero_inflation));

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n_groups; ++i) {
	if (this->groups_mask[i]) {
	    for (size_t j = 0; j < num_ecs; ++j) {
		this->log_likelihoods(groups_pos[i], j) = precalc_lls_mat(groups_pos[i], group_counts[i*num_ecs + j]);
	    }
	}
    }
  }

  void fill_ec_counts(const mSWEEP::Alignment &alignment) {
    // Fill log ec counts.
    this->log_ec_counts.resize(alignment.n_ecs(), 0);
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < alignment.n_ecs(); ++i) {
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

  LL_WOR21(const std::vector<V> &group_sizes, const mSWEEP::Alignment &alignment, const size_t n_groups, const T tol, const T frac_mu, const size_t min_hits, const T _zero_inflation) {
    this->bb_constants[0] = tol;
    this->bb_constants[1] = frac_mu;
    this->zero_inflation = _zero_inflation;
    this->from_grouped_alignment(alignment, group_sizes, n_groups, min_hits);
  }

  void from_grouped_alignment(const mSWEEP::Alignment &alignment, const std::vector<V> &group_sizes, const size_t n_groups, const size_t min_hits) {
    this->fill_ll_mat(alignment, group_sizes, n_groups, min_hits);
    this->fill_ec_counts(alignment);
  }

  void from_file(const size_t n_groups, std::istream* infile) {
    // Have to read the likelihoods into a temporary because num_ecs is not known
    std::vector<std::vector<double>> likelihoods(n_groups, std::vector<double>());
    this->log_ec_counts.resize(0);

    if (infile->good()) {
      std::string newline;
      size_t line_nr = 0;
      while (std::getline(*infile, newline)) {
	++line_nr;
	std::string part;
	std::stringstream partition(newline);
	bool ec_count_col = true;
	size_t group_id = 0;
	while (std::getline(partition, part, '\t')) {
	  if (ec_count_col) {
	    size_t ec_count = std::stol(part);
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
      for (size_t i = 0; i < n_ecs; ++i){
	size_t ec_hit_count = std::round(std::exp(log_ec_counts[i]));
	*of << ec_hit_count << '\t';
	for (size_t j = 0; j < n_groups; ++j) {
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

      size_t read_id = 1;
      for (size_t i = 0; i < n_ecs; ++i) {
	size_t ec_hit_count = std::round(std::exp(log_ec_counts[i]));
	for (size_t k = 0; k < ec_hit_count; ++k) {
	  *of << read_id << ' ';
	  *of << n_groups + 1 << ' ';
	  for (size_t j = 0; j < n_groups; ++j) {
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

  // Get the groups mask
  const std::vector<bool>& groups_considered() const override { return this->groups_mask; };
};
template <typename T>
std::unique_ptr<Likelihood<T>> ConstructAdaptiveLikelihood(const mSWEEP::Alignment &alignment, const Grouping &grouping, const T q, const T e, const size_t min_hits, const T zero_inflation) {
    size_t max_group_size = grouping.max_group_size();
    size_t n_groups = grouping.get_n_groups();
    std::unique_ptr<Likelihood<T>> log_likelihoods;
    if (max_group_size <= std::numeric_limits<uint8_t>::max()) {
	if (n_groups <= std::numeric_limits<uint8_t>::max()) {
	    log_likelihoods.reset(new mSWEEP::LL_WOR21<T, uint8_t>(static_cast<const AdaptiveGrouping<uint8_t, uint8_t>*>(&grouping)->get_sizes(), alignment,  n_groups, q, e, min_hits, zero_inflation));
	} else if (n_groups <= std::numeric_limits<uint16_t>::max()) {
	    log_likelihoods.reset(new mSWEEP::LL_WOR21<T, uint8_t>(static_cast<const AdaptiveGrouping<uint8_t, uint16_t>*>(&grouping)->get_sizes(), alignment,  n_groups, q, e, min_hits, zero_inflation));
	} else if (n_groups <= std::numeric_limits<uint32_t>::max()) {
	    log_likelihoods.reset(new mSWEEP::LL_WOR21<T, uint8_t>(static_cast<const AdaptiveGrouping<uint8_t, uint32_t>*>(&grouping)->get_sizes(), alignment,  n_groups, q, e, min_hits, zero_inflation));
	} else {
	    log_likelihoods.reset(new mSWEEP::LL_WOR21<T, uint8_t>(static_cast<const AdaptiveGrouping<uint8_t, uint64_t>*>(&grouping)->get_sizes(), alignment,  n_groups, q, e, min_hits, zero_inflation));		
	}
    } else if (max_group_size <= std::numeric_limits<uint16_t>::max()) {
	if (n_groups <= std::numeric_limits<uint8_t>::max()) {
	    log_likelihoods.reset(new mSWEEP::LL_WOR21<T, uint16_t>(static_cast<const AdaptiveGrouping<uint16_t, uint8_t>*>(&grouping)->get_sizes(), alignment,  n_groups, q, e, min_hits, zero_inflation));
	} else if (n_groups <= std::numeric_limits<uint16_t>::max()) {
	    log_likelihoods.reset(new mSWEEP::LL_WOR21<T, uint16_t>(static_cast<const AdaptiveGrouping<uint16_t, uint16_t>*>(&grouping)->get_sizes(), alignment,  n_groups, q, e, min_hits, zero_inflation));
	} else if (n_groups <= std::numeric_limits<uint32_t>::max()) {
	    log_likelihoods.reset(new mSWEEP::LL_WOR21<T, uint16_t>(static_cast<const AdaptiveGrouping<uint16_t, uint32_t>*>(&grouping)->get_sizes(), alignment,  n_groups, q, e, min_hits, zero_inflation));
	} else {
	    log_likelihoods.reset(new mSWEEP::LL_WOR21<T, uint16_t>(static_cast<const AdaptiveGrouping<uint16_t, uint64_t>*>(&grouping)->get_sizes(), alignment,  n_groups, q, e, min_hits, zero_inflation));		
	}
    } else if (max_group_size <= std::numeric_limits<uint32_t>::max()) {
	if (n_groups <= std::numeric_limits<uint8_t>::max()) {
	    log_likelihoods.reset(new mSWEEP::LL_WOR21<T, uint32_t>(static_cast<const AdaptiveGrouping<uint32_t, uint8_t>*>(&grouping)->get_sizes(), alignment,  n_groups, q, e, min_hits, zero_inflation));
	} else if (n_groups <= std::numeric_limits<uint16_t>::max()) {
	    log_likelihoods.reset(new mSWEEP::LL_WOR21<T, uint32_t>(static_cast<const AdaptiveGrouping<uint32_t, uint16_t>*>(&grouping)->get_sizes(), alignment,  n_groups, q, e, min_hits, zero_inflation));
	} else if (n_groups <= std::numeric_limits<uint32_t>::max()) {
	    log_likelihoods.reset(new mSWEEP::LL_WOR21<T, uint32_t>(static_cast<const AdaptiveGrouping<uint32_t, uint32_t>*>(&grouping)->get_sizes(), alignment,  n_groups, q, e, min_hits, zero_inflation));
	} else {
	    log_likelihoods.reset(new mSWEEP::LL_WOR21<T, uint32_t>(static_cast<const AdaptiveGrouping<uint32_t, uint64_t>*>(&grouping)->get_sizes(), alignment,  n_groups, q, e, min_hits, zero_inflation));		
	}
    } else {
	if (n_groups <= std::numeric_limits<uint8_t>::max()) {
	    log_likelihoods.reset(new mSWEEP::LL_WOR21<T, uint64_t>(static_cast<const AdaptiveGrouping<uint64_t, uint8_t>*>(&grouping)->get_sizes(), alignment,  n_groups, q, e, min_hits, zero_inflation));
	} else if (n_groups <= std::numeric_limits<uint16_t>::max()) {
	    log_likelihoods.reset(new mSWEEP::LL_WOR21<T, uint64_t>(static_cast<const AdaptiveGrouping<uint64_t, uint16_t>*>(&grouping)->get_sizes(), alignment,  n_groups, q, e, min_hits, zero_inflation));
	} else if (n_groups <= std::numeric_limits<uint32_t>::max()) {
	    log_likelihoods.reset(new mSWEEP::LL_WOR21<T, uint64_t>(static_cast<const AdaptiveGrouping<uint64_t, uint32_t>*>(&grouping)->get_sizes(), alignment,  n_groups, q, e, min_hits, zero_inflation));
	} else {
	    log_likelihoods.reset(new mSWEEP::LL_WOR21<T, uint64_t>(static_cast<const AdaptiveGrouping<uint64_t, uint64_t>*>(&grouping)->get_sizes(), alignment,  n_groups, q, e, min_hits, zero_inflation));		
	}
    }
    return log_likelihoods;
}

}

#endif

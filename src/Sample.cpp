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
#include "Sample.hpp"

#include <exception>

namespace mSWEEP {
void ConstructSample(const telescope::Alignment &alignment, const size_t bootstrap_iters, const size_t bootstrap_count, const size_t bootstrap_seed, const bool bin_reads, std::unique_ptr<Sample> &sample) {
  // Wrapper for determining which Sample type to construct.
  // Initialize Sample depending on how the alignment needs to be processed.
  if (bootstrap_iters > 0) {
    // Bootstrap mode
    bool count_provided = bootstrap_count > 0;
    if (count_provided && bin_reads) {
      sample.reset(new BinningBootstrap(alignment, bootstrap_iters, bootstrap_count, bootstrap_seed));
    } else if (count_provided && !bin_reads) {
      sample.reset(new BootstrapSample(alignment, bootstrap_iters, bootstrap_iters, bootstrap_seed));
    } else if (bin_reads) {
      sample.reset(new BinningBootstrap(alignment, bootstrap_iters, bootstrap_seed));
    } else {
      sample.reset(new BootstrapSample(alignment, bootstrap_iters, bootstrap_seed));
    }
  } else if (bin_reads) {
    sample.reset(new BinningSample(alignment));
  } else {
    sample.reset(new PlainSample(alignment));
  }
}

void Sample::count_alignments(const telescope::Alignment &alignment) {
  // Count the number of aligned reads and store in counts_total
  size_t aln_counts_total = 0;
#pragma omp parallel for schedule(static) reduction(+:aln_counts_total)
  for (size_t i = 0; i < alignment.n_ecs(); ++i) {
    aln_counts_total += alignment.reads_in_ec(i);
  }
  this->counts_total = aln_counts_total;
  this->n_reads = alignment.n_reads();
}

void Sample::write_probs(const std::vector<std::string> &cluster_indicators_to_string, std::ostream *of) {
  // Write the probability matrix to a file.
  if (of->good()) {
    *of << "ec_id" << '\t';
    size_t n_rows = this->ec_probabilities.get_rows();
    size_t n_cols = this->ec_probabilities.get_cols();
    for (size_t i = 0; i < n_rows; ++i) {
      *of << cluster_indicators_to_string[i];
      *of << (i < n_rows - 1 ? '\t' : '\n');
    }
    for (size_t i = 0; i < n_cols; ++i) {
	*of << i << '\t';
	for (size_t j = 0; j < n_rows; ++j) {
	  *of << std::exp(this->ec_probabilities(j, i));
	  *of << (j < n_rows - 1 ? '\t' : '\n');
	}
    }
    *of << std::endl;
      of->flush();
  } else {
    throw std::runtime_error("Can't write to probs file.");
  }
}

double digamma(double x) {
    double result = 0, xx, xx2, xx4;
    for ( ; x < 7; ++x)
	result -= 1/x;
    x -= 1.0/2.0;
    xx = 1.0/x;
    xx2 = xx*xx;
    xx4 = xx2*xx2;
    result += std::log(x)+(1./24.)*xx2-(7.0/960.0)*xx4+(31.0/8064.0)*xx4*xx2-(127.0/30720.0)*xx4*xx4;
    return result;
}

void Sample::dirichlet_kld(const std::vector<double> &log_ec_hit_counts) {
  size_t rows = this->get_probs().get_rows();
  size_t cols = this->get_probs().get_cols();

  std::vector<double> alphas(rows, 0.0);
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      size_t num_hits = std::round(std::exp(log_ec_hit_counts[j]));
      for (size_t k = 0; k < num_hits; ++k) {
	alphas[i] += std::exp(this->get_probs()(i, j));
      }
    }
  }

  double alpha0 = 0.0;
#pragma omp parallel for schedule(static) reduction(+:alpha0)
  for (size_t i = 0; i < rows; ++i) {
    alpha0 += alphas[i];
  }

  this->log_KLDs.resize(rows);
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < rows; ++i) {
    double log_theta = std::log(alphas[i]) - std::log(alpha0);
    double alpha_k = alphas[rows - 1];
    double alpha_j = alphas[i];
    double KLD = std::max(std::lgamma(alpha0) - std::lgamma(alpha0 - alpha_j) - std::lgamma(alpha_j) + alpha_j * (digamma(alpha_j) - digamma(alpha0)), 1e-16);
    this->log_KLDs[i] = std::log(KLD);
  }

  this->rate_run = true;
}

std::vector<double> Sample::get_rates() const {
  double max_elem = 0.0;
  // TODO pragma with custom reduction to find maximum
  for (size_t i = 0; i < this->log_KLDs.size(); ++i) {
    max_elem = (max_elem > this->log_KLDs[i] ? max_elem : this->log_KLDs[i]);
  }
  double tmp_sum = 0.0;
#pragma omp parallel for schedule(static) reduction(+:tmp_sum)
  for (size_t i = 0; i < this->log_KLDs.size(); ++i) {
    tmp_sum += std::exp(this->log_KLDs[i] - max_elem);
  }
  double log_KLDs_sum = std::log(tmp_sum) + max_elem;

  std::vector<double> RATE(this->log_KLDs.size());
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < this->log_KLDs.size(); ++i) {
    RATE[i] = std::exp(this->log_KLDs[i] - log_KLDs_sum);
  }
  return RATE;
}

void Sample::write_probs2(const std::vector<std::string> &estimated_indicators_to_string,
			  const std::vector<std::string> &zero_indicators_to_string, std::ostream *of) {
  // Write the probability matrix to a file.
  if (of->good()) {
    *of << "ec_id" << '\t';
    size_t n_rows = estimated_indicators_to_string.size() + zero_indicators_to_string.size();
    size_t n_cols = this->ec_probabilities.get_cols();
    for (size_t i = 0; i < n_rows; ++i) {
	if (i < estimated_indicators_to_string.size()) {
	    *of << estimated_indicators_to_string[i];
	    *of << (i < n_rows - 1 ? '\t' : '\n');
	} else {
	    *of << zero_indicators_to_string[i - estimated_indicators_to_string.size()];
	    *of << (i < n_rows - 1 ? '\t' : '\n');
	}
    }
    for (size_t i = 0; i < n_cols; ++i) {
	*of << i << '\t';
	for (size_t j = 0; j < n_rows; ++j) {
	    if (j < estimated_indicators_to_string.size()) {
		*of << std::exp(this->ec_probabilities(j, i));
	    } else {
		*of << (double)0.0;
	    }
	  *of << (j < n_rows - 1 ? '\t' : '\n');
	}
    }
    *of << std::endl;
      of->flush();
  } else {
    throw std::runtime_error("Can't write to probs file.");
  }
}
}

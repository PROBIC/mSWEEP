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
#ifndef MSWEEP_SAMPLE_HPP
#define MSWEEP_SAMPLE_HPP

#include <cstddef>
#include <vector>
#include <string>
#include <fstream>
#include <random>
#include <memory>

#include "Matrix.hpp"
#include "telescope.hpp"

namespace mSWEEP {
class Sample {
private:
  size_t n_reads;
  size_t counts_total;
  seamat::DenseMatrix<double> ec_probabilities;

  // Relative abundance estimate scoring via RATEs
  bool rate_run = false;
  std::vector<double> log_KLDs;

protected:
  void count_alignments(const telescope::Alignment &alignment);

public:
  // Virtual functions
  // Store the relative abundances for later.
  virtual void store_abundances(const std::vector<double> &abundances) =0;

  // Getters
  virtual const std::vector<double>& get_abundances() const =0;
  // Write the relative abundances
  virtual void write_abundances(const std::vector<std::string> &group_names, std::ostream *of) const =0;
  virtual void write_abundances2(const std::vector<std::string> &estimated_group_names,
				 const std::vector<std::string> &zero_group_names, std::ostream *of) const =0;

  // Non-virtuals
  // Store equivalence class probabilities
  void store_probs(const seamat::DenseMatrix<double> &probs) { this->ec_probabilities = std::move(probs); }

  // Write equivalence class probabilities
  void write_probs(const std::vector<std::string> &cluster_indicators_to_string, std::ostream *of);
  void write_probs2(const std::vector<std::string> &cluster_indicators_to_string,
		    const std::vector<std::string> &zero_indicators_to_string, std::ostream *of);

  // Calculate KLDs from probs
  void dirichlet_kld(const std::vector<double> &log_ec_hit_counts);

  // Getters
  size_t get_counts_total() const { return this->counts_total; };
  size_t get_n_reads() const { return this->n_reads; };

  const seamat::DenseMatrix<double>& get_probs() const { return this->ec_probabilities; }
  const std::vector<double>& get_log_klds() const { return this->log_KLDs; }
  std::vector<double> get_rates() const;

  size_t get_n_ecs() const { return this->ec_probabilities.get_rows(); }
  size_t get_n_refs() const { return this->ec_probabilities.get_cols(); }
  size_t get_rate_run() const { return this->rate_run; }
};

class Binning {
private:
  // Need to store the read assignments to equivalence classes if also binning
  std::vector<std::vector<uint32_t>> aligned_reads;

protected:
  // This class should never be initialized.
  Binning() = default;

  // Function for derived classes to set aligned_reads.
  void store_aligned_reads(const std::vector<std::vector<uint32_t>> &_aligned_reads) { this->aligned_reads = _aligned_reads; }

public:
  // Getters
  const std::vector<std::vector<uint32_t>>& get_aligned_reads() const { return this->aligned_reads; }

};

class PlainSample : public Sample {
private:
  std::vector<double> relative_abundances;

public:
  PlainSample() = default;

  PlainSample(const telescope::Alignment &alignment) {
    this->count_alignments(alignment);
  }

  // Store relative abundances for writing
  void store_abundances(const std::vector<double> &abundances) override { this->relative_abundances = std::move(abundances); }

  // Write the relative abundances
  void write_abundances(const std::vector<std::string> &group_names, std::ostream *of) const override;
  void write_abundances2(const std::vector<std::string> &estimated_group_names,
			 const std::vector<std::string> &zero_group_names, std::ostream *of) const override;

  // Getters
  const std::vector<double>& get_abundances() const override { return this->relative_abundances; }

};

class BinningSample : public PlainSample, public Binning {
public:
  BinningSample() = default;

  BinningSample(const telescope::Alignment &alignment) {
    this->count_alignments(alignment);
    this->store_aligned_reads(alignment.get_aligned_reads());
  }

};

class BootstrapSample : public Sample {
private:
  std::mt19937_64 gen;
  std::discrete_distribution<uint32_t> ec_distribution;

  // Max number of bootstrap iterations
  size_t iters;

  // Number of ecs to bootstrap
  size_t bootstrap_count;

  // Need to store this for resampling counts.
  size_t num_ecs;

  // First element has the relative abundances without bootstrapping.
  std::vector<std::vector<double>> bootstrap_results;

  // Set all variables required to bootstrap the ec_counts later
  void init_bootstrap(const telescope::Alignment &alignment);

protected:
  void construct(const telescope::Alignment &alignment, const size_t _iters, const int32_t seed, const size_t bootstrap_count=0);

public:
  BootstrapSample() = default;

  // Set seed in constructor
  BootstrapSample(const telescope::Alignment &alignment, const size_t _iters, const int32_t seed) {
    this->construct(alignment, _iters, seed);
  }
  BootstrapSample(const telescope::Alignment &alignment, const size_t _iters, const size_t _bootstrap_count, const int32_t seed) {
    this->construct(alignment, _iters, seed, _bootstrap_count);
  }

  // Resample the equivalence class counts
  std::vector<double> resample_counts();

  // Store relative abundances in bootstrap_results
  void store_abundances(const std::vector<double> &abundances) override { this->bootstrap_results.emplace_back(std::move(abundances)); }

  // Write the bootstrap results
  void write_abundances(const std::vector<std::string> &group_names, std::ostream *os) const override;
  void write_abundances2(const std::vector<std::string> &estimated_group_names,
			 const std::vector<std::string> &zero_group_names, std::ostream *of) const override;

  // Getters
  const std::vector<double>& get_abundances() const override { return this->bootstrap_results[0]; }

};

class BinningBootstrap : public BootstrapSample, public Binning {
public:
  BinningBootstrap(const telescope::Alignment &alignment, const size_t _iters, const int32_t seed) {
    this->construct(alignment, _iters, seed);
    this->store_aligned_reads(alignment.get_aligned_reads());
  }
  BinningBootstrap(const telescope::Alignment &alignment, const size_t _iters, const size_t _bootstrap_count, const int32_t seed) {
    this->construct(alignment, _iters, seed, _bootstrap_count);
    this->store_aligned_reads(alignment.get_aligned_reads());
  }

};

void ConstructSample(const telescope::Alignment &alignment, const size_t bootstrap_iters, const size_t bootstrap_count, const size_t bootstrap_seed, const bool bin_reads, std::unique_ptr<Sample> &sample);

}

#endif

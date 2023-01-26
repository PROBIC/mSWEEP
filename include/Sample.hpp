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
#include <random>

#include "telescope.hpp"

class Sample {
private:
  uint32_t counts_total;

protected:
  void count_alignments(const telescope::Alignment &alignment) {
    // Count the number of aligned reads and store in counts_total
    uint32_t aln_counts_total = 0;
#pragma omp parallel for schedule(static) reduction(+:aln_counts_total)
    for (uint32_t i = 0; i < alignment.n_ecs(); ++i) {
      aln_counts_total += alignment.reads_in_ec(i);
    }
    this->counts_total = aln_counts_total;
  }

public:

  // Virtual functions
  // Store the relative abundances for later.
  virtual void store_abundances(const std::vector<double> &abundances) =0;
  // Getters
  virtual const std::vector<double>& get_abundances() const =0;
  // Write the relative abundances
  virtual void write_abundances(const std::vector<std::string> &group_names, std::ostream *of) const =0;

  // Non-virtuals
  // Getters
  uint32_t get_counts_total() const { return this->counts_total; };

};

class PlainSample : public Sample {
private:
  std::vector<double> relative_abundances;

public:
  PlainSample() = default;
  PlainSample(const telescope::GroupedAlignment &alignment) {
    this->count_alignments(alignment);
  }

  // Store relative abundances for writing
  void store_abundances(const std::vector<double> &abundances) override { this->relative_abundances = std::move(abundances); }

  // Write the relative abundances
  void write_abundances(const std::vector<std::string> &group_names, std::ostream *of) const override;

  // Getters
  const std::vector<double>& get_abundances() const override { return this->relative_abundances; }

};

class BinningSample : public PlainSample {
private:
  // Need to store the read assignments to equivalence classes if also binning
  std::vector<std::vector<uint32_t>> aligned_reads;

public:
  BinningSample(const telescope::GroupedAlignment &alignment) {
    this->count_alignments(alignment);
    this->aligned_reads = alignment.get_aligned_reads();
  }

  // Getters
  const std::vector<std::vector<uint32_t>>& get_aligned_reads() const { return this->aligned_reads; }

};

class BootstrapSample : public Sample {
private:
  std::mt19937_64 gen;
  std::discrete_distribution<uint32_t> ec_distribution;

  // Max number of bootstrap iterations
  size_t iters;

  // Need to store this for resampling counts.
  size_t num_ecs;

  // First element has the relative abundances without bootstrapping.
  std::vector<std::vector<double>> bootstrap_results;

  // Set all variables required to bootstrap the ec_counts later
  void init_bootstrap(const telescope::GroupedAlignment &alignment);

public:
  // Set seed in constructor
  BootstrapSample(const telescope::GroupedAlignment &alignment, const size_t _iters, const int32_t seed);

  // Resample the equivalence class counts
  std::vector<double> resample_counts(const uint32_t how_many);

  // Store relative abundances in bootstrap_results
  void store_abundances(const std::vector<double> &abundances) override { this->bootstrap_results.emplace_back(std::move(abundances)); }

  // Write the bootstrap results
  void write_abundances(const std::vector<std::string> &group_names, std::ostream *os) const override;

  // Getters
  const std::vector<double>& get_abundances() const override { return this->bootstrap_results[0]; }

};

#endif

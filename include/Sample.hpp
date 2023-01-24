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

  // Need to store the read assignments to equivalence classes if also binning
  std::vector<std::vector<uint32_t>> aligned_reads;

public:
  Sample() = default;
  Sample(const telescope::GroupedAlignment &alignment, bool bin_reads) {
    uint32_t aln_counts_total = 0;
#pragma omp parallel for schedule(static) reduction(+:aln_counts_total)
    for (uint32_t i = 0; i < alignment.n_ecs(); ++i) {
      aln_counts_total += alignment.reads_in_ec(i);
    }
    counts_total = aln_counts_total;

    if (bin_reads) {
      this->aligned_reads = alignment.get_aligned_reads();
    }
  }

  // Getters
  uint32_t get_counts_total() const { return this->counts_total; };
  const std::vector<std::vector<uint32_t>>& get_aligned_reads() const { return this->aligned_reads; }

};

class BootstrapSample : public Sample {
private:
  std::mt19937_64 gen;
  std::discrete_distribution<uint32_t> ec_distribution;

  // Need to store this for resampling counts.
  size_t num_ecs;

  // First element has the relative abundances without bootstrapping.
  std::vector<std::vector<double>> bootstrap_results;

  // Set all variables required to bootstrap the ec_counts later
  void init_bootstrap(const telescope::GroupedAlignment &alignment);

public:
  // Set seed in constructor
  BootstrapSample(const telescope::GroupedAlignment &alignment, const int32_t seed);

  // Resample the equivalence class counts
  std::vector<double> resample_counts(const uint32_t how_many);

  // Store relative abundances in bootstrap_results
  void move_abundances(const std::vector<double> &relative_abundances) { this->bootstrap_results.emplace_back(std::move(relative_abundances)); }

  // Get the results
  const std::vector<std::vector<double>>& get_results() const { return this->bootstrap_results; }

};

#endif

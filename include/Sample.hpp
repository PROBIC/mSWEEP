#ifndef MSWEEP_SAMPLE_HPP
#define MSWEEP_SAMPLE_HPP

#include <cstddef>
#include <string>
#include <vector>
#include <fstream>
#include <random>

#include "telescope.hpp"
#include "Matrix.hpp"

#include "Grouping.hpp"
#include "Reference.hpp"

class Sample {
private:
  uint32_t counts_total;

public:
  Sample() = default;
  Sample(const telescope::GroupedAlignment &alignment) {
    uint32_t aln_counts_total = 0;
#pragma omp parallel for schedule(static) reduction(+:aln_counts_total)
    for (uint32_t i = 0; i < alignment.n_ecs(); ++i) {
      aln_counts_total += alignment.reads_in_ec(i);
    }
    counts_total = aln_counts_total;
  }

  telescope::GroupedAlignment pseudos;

  // Getters
  uint32_t get_counts_total() const { return this->counts_total; };
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

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
#include "parse_arguments.hpp"

class Sample {
private:
  uint32_t m_num_ecs;
  std::string cell_id;

public:
  Sample() = default;
  Sample(const Reference &reference)
    : pseudos(telescope::GroupedAlignment(reference.get_n_refs(), reference.get_grouping(0).get_n_groups(), reference.get_group_indicators(0)))
  { };

  uint32_t counts_total;

  seamat::IndexMatrix<double, uint16_t> ll_mat;
  seamat::DenseMatrix<double> ec_probs;
  std::vector<double> log_ec_counts;
  std::vector<double> relative_abundances;

  // Alignments class from telescope
  telescope::GroupedAlignment pseudos;

  // Calculate log_ec_counts and counts_total.
  void process_aln(const bool bootstrap_mode);

  // Count the number of pseudoalignments in groups defined by the given indicators.
  std::vector<uint16_t> group_counts(const std::vector<uint32_t> indicators, const uint32_t ec_id, const uint32_t n_groups) const;

  // Write estimated relative abundances
  void write_abundances(const std::vector<std::string> &cluster_indicators_to_string, std::ostream &of) const;
  // Write estimated read-reference posterior probabilities (gamma_Z)
  void write_probabilities(const std::vector<std::string> &cluster_indicators_to_string, std::ostream &outfile) const;
  // Write likelihoods
  void write_likelihood(const uint32_t n_groups, std::ostream &of) const;
  // Write likelihoods in BitSeq-compatible format
  void write_likelihood_bitseq(const uint32_t n_groups, std::ostream &of) const;

  // Read in the likelihoods from a file
  void read_likelihood(const Grouping &grouping, std::istream &infile);

  // Getters
  std::string cell_name() const { return cell_id; };
  uint32_t num_ecs() const { return m_num_ecs; };
  uint32_t get_counts_total() const { return this->counts_total; };
};

class BootstrapSample : public Sample {
private:
  std::mt19937_64 gen;
  std::discrete_distribution<uint32_t> ec_distribution;

  // Run estimation and add results to relative_abundances
  void bootstrap_iter(const std::vector<double> &resampled_log_ec_counts,
		      const std::vector<double> &alpha0, const double tolerance,
		      const uint16_t max_iters);

public:
  // Set seed in constructor
  BootstrapSample(const int32_t seed);

  std::vector<std::vector<double>> bootstrap_results;

  // Resample the equivalence class counts
  std::vector<double> resample_counts(const uint32_t how_many);

  void init_bootstrap();

  // Estimate the mixture components with bootstrap iterations
  void estimate_abundances(const Arguments &args);

  void write_bootstrap(const std::vector<std::string> &cluster_indicators_to_string,
		       const uint16_t iters, std::ostream &of) const;
  void bootstrap_ec_counts(const Arguments &args);

};

#endif

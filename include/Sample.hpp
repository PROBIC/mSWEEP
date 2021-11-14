#ifndef MSWEEP_SAMPLE_HPP
#define MSWEEP_SAMPLE_HPP

#include <string>
#include <vector>
#include <fstream>
#include <random>

#include "telescope.hpp"

#include "Matrix.hpp"
#include "Reference.hpp"
#include "parse_arguments.hpp"

class Sample {
private:
  uint32_t m_num_ecs;
  std::string cell_id;

protected:
  uint32_t counts_total;

public:
  rcgpar::Matrix<double> ec_probs;
  rcgpar::Matrix<double> ll_mat;
  std::vector<double> log_ec_counts;

  // Count the number of pseudoalignments in groups defined by the given indicators.
  std::vector<uint16_t> group_counts(const std::vector<uint32_t> indicators, const uint32_t ec_id, const uint32_t n_groups) const;

  KallistoAlignment pseudos;

  // Calculate log_ec_counts and counts_total.
  void process_aln(const bool bootstrap_mode);

  // Write estimated relative abundances
  void write_abundances(const std::vector<std::string> &cluster_indicators_to_string, std::string outfile) const;
  // Write estimated read-reference posterior probabilities (gamma_Z)
  void write_probabilities(const std::vector<std::string> &cluster_indicators_to_string, std::ostream &outfile) const;
  // Write likelihoods
  void write_likelihood(const bool gzip_output, const uint32_t n_groups, std::string outfile) const;
  void write_likelihood_bitseq(const bool gzip_output, const uint32_t n_groups, std::string outfile) const;
  // Getters
  std::string cell_name() const { return cell_id; };
  uint32_t num_ecs() const { return m_num_ecs; };

  // Read in the likelihoods from a file
  void read_likelihood(const Grouping &grouping, std::istream &infile);
};

class BootstrapSample : public Sample {
private:
  std::mt19937_64 gen;
  std::discrete_distribution<uint32_t> ec_distribution;
  std::vector<std::vector<double>> relative_abundances;

  // Run estimation and add results to relative_abundances
  void bootstrap_iter(const std::vector<double> &resampled_log_ec_counts,
		      const std::vector<double> &alpha0, const double tolerance,
		      const uint16_t max_iters);

  // Resample the equivalence class counts
  std::vector<double> resample_counts(const uint32_t how_many, std::mt19937_64 &rng);

public:
  // Set seed in constructor
  BootstrapSample(const int32_t seed);

  void write_bootstrap(const std::vector<std::string> &cluster_indicators_to_string, std::string outfile, const unsigned iters, const bool batch_mode) const;
  void bootstrap_abundances(const Grouping &grouping, const Arguments &args);

};

#endif

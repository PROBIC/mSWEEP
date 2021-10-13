#ifndef MSWEEP_SAMPLE_HPP
#define MSWEEP_SAMPLE_HPP

#include <string>
#include <vector>
#include <fstream>
#include <random>

#include "telescope.hpp"

#include "matrix.hpp"
#include "Reference.hpp"
#include "parse_arguments.hpp"

class VSample {
public:
  virtual void read_themisto(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &strands) =0;
  virtual void read_kallisto(const uint32_t n_refs, std::istream &tsv_file, std::istream &ec_file) =0;
};

class Sample : public VSample{
private:
  uint32_t m_num_refs;
  uint32_t m_num_ecs;
  std::string cell_id;

  // Count the number of pseudoalignments in groups defined by the given indicators.
  std::vector<uint16_t> group_counts(const std::vector<uint32_t> indicators, const uint32_t ec_id, const uint32_t n_groups) const;

  // Free the memory taken by ec_configs
  void clear_configs() { pseudos.ec_configs.clear(); }

protected:
  // Calculate log_ec_counts and counts_total.
  void process_aln();

  KallistoAlignment pseudos;
  uint32_t counts_total;

public:
  Matrix<double> ec_probs;
  Matrix<double> ll_mat;
  std::vector<std::vector<uint16_t>> counts;
  std::vector<double> log_ec_counts;

  // Retrieve relative abundances from the ec_probs matrix.
  std::vector<double> group_abundances() const;
  // Write estimated relative abundances
  void write_abundances(const std::vector<std::string> &cluster_indicators_to_string, std::string outfile) const;
  // Write estimated read-reference posterior probabilities (gamma_Z)
  void write_probabilities(const std::vector<std::string> &cluster_indicators_to_string, const bool gzip_probs, std::ostream &outfile) const;
  // Write likelihoods
  void write_likelihood(const bool gzip_output, const uint32_t n_groups, std::string outfile) const;
  void write_likelihood_bitseq(const bool gzip_output, const uint32_t n_groups, std::string outfile) const;
  // Getters
  std::string cell_name() const { return cell_id; };
  uint32_t num_ecs() const { return m_num_ecs; };
  uint32_t total_counts() const { return counts_total; };

  // Read Themisto or kallisto pseudoalignments
  void read_themisto(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &strands) override;
  void read_kallisto(const uint32_t n_refs, std::istream &tsv_file, std::istream &ec_file) override;
  // Fill the likelihood matrix
  void CalcLikelihood(const Grouping &grouping, const std::vector<uint32_t> &group_indicators);
};

class BootstrapSample : public Sample {
private:
  std::discrete_distribution<uint32_t> ec_distribution;
  std::vector<std::vector<double>> relative_abundances;

  // Run estimation and add results to relative_abundances
  void BootstrapIter(const std::vector<double> &alpha0, const double tolerance, const uint16_t max_iters);
  // Initialize ec_distributino and ll_mat for bootstrapping
  void InitBootstrap(const Grouping &grouping, const std::vector<uint32_t> &group_indicators);
  // Resample the equivalence class counts
  void ResampleCounts(const uint32_t how_many, std::mt19937_64 &rng);

public:
  void WriteBootstrap(const std::vector<std::string> &cluster_indicators_to_string, std::string &outfile, const unsigned iters, const bool batch_mode) const;
  void BootstrapAbundances(const Reference &reference, const Arguments &args);

  // Read in pseudoalignments but do not free the memory used by storing the equivalence class counts.
  void read_themisto(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &strands) override;
  void read_kallisto(const uint32_t n_refs, std::istream &tsv_file, std::istream &ec_file) override;
};

#endif

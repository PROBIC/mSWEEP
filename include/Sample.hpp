#ifndef MSWEEP_SAMPLE_HPP
#define MSWEEP_SAMPLE_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <random>

#include "telescope.hpp"

#include "matrix.hpp"
#include "Reference.hpp"
#include "parse_arguments.hpp"

class Sample {
private:
  uint32_t m_num_refs;

  // Calculate log_ec_counts and counts_total.
  void process_aln();
protected:
  std::string cell_id;
  // Note the following two should be member variables in BootstrapSample, todo.
  std::vector<std::vector<double>> relative_abundances;
  std::discrete_distribution<uint32_t> ec_distribution;

  KallistoAlignment pseudos;
  uint32_t counts_total;
  uint32_t m_num_ecs;
public:
  Matrix<double> ll_mat;
  Matrix<double> ec_probs;
  std::vector<std::vector<uint16_t>> counts;
  std::vector<double> log_ec_counts;

  // Read Themisto or kallisto pseudoalignments
  void read_themisto(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &strands);
  void read_kallisto(const uint32_t n_refs, std::istream &tsv_file, std::istream &ec_file);

  // Retrieve relative abundances from the ec_probs matrix.
  std::vector<double> group_abundances() const;

  // Write estimated relative abundances
  void write_abundances(const std::vector<std::string> &cluster_indicators_to_string, std::string outfile) const;
  // Write estimated read-reference posterior probabilities (gamma_Z)
  void write_probabilities(const std::vector<std::string> &cluster_indicators_to_string, const bool gzip_probs, std::ostream &outfile) const;

  // Count the number of pseudoalignments in groups defined by the given indicators.
  std::vector<uint16_t> group_counts(const std::vector<uint32_t> indicators, const uint32_t ec_id, const uint32_t n_groups) const;

  void clear_configs() { pseudos.ec_configs.clear(); }
  // Getters
  uint32_t num_ecs() const { return m_num_ecs; };
  const std::string &cell_name() const { return cell_id; };
  const uint32_t &total_counts() const { return counts_total; };
};

class BootstrapSample : public Sample {
private:
  // Run estimation and add results to relative_abundances
  void BootstrapIter(const std::vector<double> &alpha0, const double tolerance, const uint16_t max_iters);
  // Initialize ec_distributino and ll_mat for bootstrapping
  void InitBootstrap(Grouping &grouping);
  // Resample the equivalence class counts
  void ResampleCounts(std::mt19937_64 &rng);
public:
  void BootstrapAbundances(Reference &reference, Arguments &args);
  void WriteBootstrap(const std::vector<std::string> &cluster_indicators_to_string, std::string &outfile, const unsigned iters, const bool batch_mode);
};

#endif

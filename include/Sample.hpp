#ifndef MSWEEP_SAMPLE_HPP
#define MSWEEP_SAMPLE_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <random>

#include "telescope/include/telescope.hpp"

#include "matrix.hpp"
#include "Reference.hpp"

class Sample {
protected:
  KallistoAlignment aln;
public:
  unsigned m_num_ecs;
  short unsigned m_num_refs;
  std::string cell_id;
  //  std::vector<std::vector<bool>> ec_configs;
  //  std::vector<long unsigned> ec_ids;
  //  std::vector<unsigned> ec_counts;
  std::vector<double> log_ec_counts;
  unsigned counts_total;
  Matrix<double> ec_probs;
  std::vector<std::vector<short unsigned>> counts;

  Sample() = default;
  ~Sample() = default;
  Sample(const Sample &sample) = default;
  Sample& operator=(const Sample &t) = default;

  // Parse the ec_configs and ec_counts after they've been filled.
  void process_aln(const bool create_ids);

  // Retrieve relative abundances from the ec_probs matrix.
  std::vector<double> group_abundances() const;

  // Writer functions for the contents.
  void write_probabilities(const std::vector<std::string> &cluster_indicators_to_string, const bool gzip_probs, std::ostream &outfile) const;
  void write_abundances(const std::vector<std::string> &cluster_indicators_to_string, std::string outfile) const;

  // Count the number of pseudoalignments in groups defined by the given indicators.
  std::vector<short unsigned> group_counts(const std::vector<unsigned short> indicators, const unsigned ec_id, const uint16_t n_groups) const;

  void clear_configs() { access_aln()->access_aln()->ec_configs.clear(); }
  void clear_counts() { access_aln()->access_aln()->ec_counts.clear(); }
  void clear_ids() { aln.clear_ids(); }

  // Getters
  unsigned num_ecs() const { return m_num_ecs; };
  const std::string &cell_name() const { return cell_id; };
  const unsigned &total_counts() const { return counts_total; };
  const std::vector<std::vector<bool>> &ec_configs() const { return aln.get_ec_configs(); }
  const std::vector<uint32_t> &ec_ids() const { return aln.get_ec_ids(); }
  const std::vector<uint32_t> &ec_counts() const { return aln.get_ec_counts(); }

  KallistoAlignment* access_aln() { return &aln; }
};

class SampleBS : public Sample {
private:
  std::discrete_distribution<unsigned> ec_distribution;
public:
  Matrix<double> ll_mat = Matrix<double>(0, 0, 0.0);
  void init_bootstrap(Grouping &grouping);
  void resample_counts(std::mt19937_64 &rng);
  void write_bootstrap(const std::vector<std::string> &cluster_indicators_to_string, std::string outfile, unsigned iters);
  std::unordered_map<unsigned, std::vector<double>> bootstrap_abundances;
};

#endif

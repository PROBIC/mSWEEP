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
private:

  // Bootstrapping variables
  std::discrete_distribution<long unsigned> ec_distribution;

public:
  short unsigned n_groups;
  unsigned m_num_ecs;
  short unsigned m_num_refs;
  std::string cell_id;
  std::vector<std::vector<bool>> ec_configs;
  std::vector<long unsigned> ec_ids;
  std::vector<long unsigned> ec_counts;
  long unsigned counts_total;
  Matrix<double> ec_probs = Matrix<double>(0, 0, 0.0);
  // Optional storage for likelihood, used in bootstrap
  Matrix<double> ll_mat = Matrix<double>(0, 0, 0.0);
  std::vector<std::vector<short unsigned>> group_hitcounts;

  // Bootstrap results
  std::unordered_map<unsigned, std::vector<double>> bootstrap_abundances;

  Sample() = default;
  Sample(std::string cell_id_p, std::vector<long unsigned> ec_ids_p, std::vector<long unsigned> ec_counts_p, long unsigned counts_total_p, std::vector<std::vector<bool>> ec_configs_p);
  Sample(KAlignment converted_aln);

  // Retrieve relative abundances from the ec_probs matrix.
  std::vector<double> group_abundances() const;

  // Writer functions for the contents.
  void write_probabilities(const std::vector<std::string> &cluster_indicators_to_string, const bool gzip_probs, std::ostream &outfile) const;
  void write_abundances(const std::vector<std::string> &cluster_indicators_to_string, std::string outfile) const;
  void write_bootstrap(const std::vector<std::string> &cluster_indicators_to_string, std::string outfile, unsigned iters);

  // Count the number of pseudoalignments in groups defined by the given indicators.
  std::vector<short unsigned> group_counts(const std::vector<unsigned short> indicators, unsigned ec_id) const;
  // Return the number of pseudoalignments in a given group
  std::vector<bool> group_counts(unsigned ec_id_pos) const { return this->ec_configs[ec_id_pos]; };

  short unsigned group_counts(unsigned ec_id_pos, unsigned group_id) const { return this->ec_configs[ec_id_pos][group_id]; }

  // Initialize bootstrapping variables
  void init_bootstrap(Grouping &grouping);
  // Resample the pseudoalignment counts
  void resample_counts(std::mt19937_64 &rng);

  // Getters
  unsigned num_ecs() const { return this->ec_ids.size(); };
  const std::string &cell_name() const { return this->cell_id; };
  const long unsigned &total_counts() const { return this->counts_total; };
};

struct BootstrapResults {
  std::unordered_map<std::string, std::pair<unsigned, std::vector<std::vector<double>>>> results;

  //  void at(std::string key) { return this->results.at(key); };
  void insert(std::string key, unsigned counts, std::vector<std::vector<double>> abundances) { this->results.insert(std::make_pair(key, std::make_pair(counts, abundances))); };
  void insert_iter(std::string key, std::vector<double> iter) { this->results.at(key).second.emplace_back(iter); };
  std::unordered_map<std::string, std::pair<unsigned, std::vector<std::vector<double>>>> get() const { return this->results; };
};

#endif

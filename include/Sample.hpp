#ifndef SAMPLE_H
#define SAMPLE_H

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <fstream>
#include <iostream>
#include <random>

#include "matrix.hpp"

class Sample {
private:
  std::vector<long unsigned> ec_ids;
  std::shared_ptr<std::unordered_map<long unsigned, std::vector<bool>>> ec_configs;
  std::string cell_id;
  long unsigned counts_total;

  // Bootstrapping variables
  std::discrete_distribution<long unsigned> ec_distribution;
  
public:
  std::vector<long unsigned> ec_counts;
  Matrix<double> ec_probs = Matrix<double>(0, 0, 0.0);

  Sample(std::string cell_id_p, std::vector<long unsigned> ec_ids_p, std::vector<long unsigned> ec_counts_p, long unsigned counts_total_p, std::shared_ptr<std::unordered_map<long unsigned, std::vector<bool>>> ec_configs_p);

  // Retrieve relative abundances from the ec_probs matrix.
  std::vector<double> group_abundances() const;

  // Writer functions for the contents.
  void write_probabilities(const std::vector<std::string> &cluster_indicators_to_string, std::string outfile) const;
  void write_abundances(const std::vector<std::string> &cluster_indicators_to_string, std::string outfile) const;

  // Count the number of pseudoalignments in groups defined by the given indicators.
  std::vector<unsigned> group_counts(const std::vector<signed> &indicators, unsigned n_groups, unsigned ec_id_pos) const;

  // Initialize bootstrapping variables
  void init_bootstrap() { this->ec_distribution = std::discrete_distribution<long unsigned>(this->ec_counts.begin(), this->ec_counts.end()); };
  // Resample the pseudoalignment counts
  void resample_counts(std::mt19937_64 &rng);

  // Getters
  unsigned num_ecs() const { return this->ec_ids.size(); };
  const std::string &cell_name() const { return this->cell_id; };
  const long unsigned &total_counts() const { return this->counts_total; };
};

#endif

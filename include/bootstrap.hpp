#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

#include <vector>
#include <string>
#include <unordered_map>

#include "thread_pool.hpp"
#include "Sample.hpp"
#include "parse_arguments.hpp"

struct BootstrapResults {
  std::unordered_map<std::string, std::pair<unsigned, std::vector<std::vector<double>>>> results;

  void insert(std::string key, unsigned counts, std::vector<std::vector<double>> abundances) { this->results.insert(std::make_pair(key, std::make_pair(counts, abundances))); };
  void insert_iter(std::string key, std::vector<double> iter) { this->results.at(key).second.emplace_back(iter); };
  std::unordered_map<std::string, std::pair<unsigned, std::vector<std::vector<double>>>> get() const { return this->results; };
};

void write_bootstrap(const std::vector<std::string> &cluster_indicators_to_string, const std::vector<std::vector<double>> &abundances, std::string &outfile, unsigned iters, unsigned counts_total);
BootstrapResults bootstrap_abundances(const std::vector<Sample> &bitfields, Reference &reference, ThreadPool &pool, Arguments &args);

#endif

#include "Sample.hpp"

#include <memory>

#include "likelihood.hpp"
#include "version.h"

Sample::Sample(std::string cell_id_p, std::vector<long unsigned> ec_ids_p, std::vector<long unsigned> ec_counts_p, long unsigned counts_total_p, std::shared_ptr<std::unordered_map<long unsigned, std::vector<bool>>> ec_configs_p) {
  this->cell_id = cell_id_p;
  this->ec_ids = ec_ids_p;
  this->ec_counts = ec_counts_p;
  this->counts_total = counts_total_p;
  this->ec_configs = ec_configs_p;
}

Sample::Sample(KAlignment converted_aln) {
  this->cell_id = "";
  this->ec_ids.resize(converted_aln.ecs.size());
  this->ec_counts.resize(converted_aln.ecs.size());
  this->ec_configs.reset(new std::unordered_map<long unsigned, std::vector<bool>>);
  this->ec_configs->reserve(converted_aln.ecs.bucket_count());
  this->counts_total = 0;
  size_t i = 0;
  for (auto kv : converted_aln.ecs) {
    this->ec_ids[i] = i;
    this->ec_counts[i] = kv.second.count;
    this->counts_total += kv.second.count;
    this->ec_configs->insert(std::make_pair(i, kv.first));
    ++i;
  }
}

std::vector<double> Sample::group_abundances() const {
  // Calculate the relative abundances of the
  // reference groups from the ec_probs matrix
  std::vector<double> thetas(this->ec_probs.get_rows(), 0.0);
  for (unsigned i = 0; i < this->ec_probs.get_rows(); ++i) {
    for (unsigned j = 0; j < this->ec_probs.get_cols(); ++j) {
      thetas[i] += this->ec_probs(i, j) * this->ec_counts[j];
    }
    thetas[i] /= this->counts_total;
  }
  return thetas;
}

std::vector<unsigned> Sample::group_counts(const std::vector<signed> &indicators, unsigned n_groups, unsigned ec_id_pos) const {
  long unsigned ec_id = this->ec_ids[ec_id_pos];
  std::vector<bool> bitset = (*this->ec_configs)[ec_id];

  std::vector<unsigned> read_hitcounts(n_groups);
  for (size_t j = 0; j < bitset.size(); ++j) {
    read_hitcounts[indicators[j]] += bitset[j];
  }
  return read_hitcounts;
}

void Sample::init_bootstrap(Grouping &grouping) {
  this->ec_distribution = std::discrete_distribution<long unsigned>(this->ec_counts.begin(), this->ec_counts.end());
  this->ll_mat = likelihood_array_mat(*this, grouping);
}

void Sample::resample_counts(std::mt19937_64 &generator) {
  std::vector<long unsigned> new_counts(this->ec_counts.size());
  for (long unsigned i = 0; i < this->counts_total; ++i) {
    long unsigned ec_id = this->ec_distribution(generator);
    new_counts[ec_id] += 1;
  }
  this->ec_counts = new_counts;
}

void Sample::write_probabilities(const std::vector<std::string> &cluster_indicators_to_string, const bool gzip_probs, std::ostream &of) const {
  // Write the probability matrix to a file.
  if (of.good()) {
    of << "ec_id" << ',';
    for (unsigned i = 0; i < this->ec_probs.get_rows(); ++i) {
      of << cluster_indicators_to_string[i];
      of << (i < this->ec_probs.get_rows() - 1 ? ',' : '\n');
    }
    for (unsigned i = 0; i < this->ec_probs.get_cols(); ++i) {
      of << this->ec_ids[i] << ',';
      for (unsigned j = 0; j < this->ec_probs.get_rows(); ++j) {
	of << this->ec_probs(j, i);
	of << (j < this->ec_probs.get_rows() - 1 ? ',' : '\n');
      }
    }
  }
  of << std::endl;
  of.flush();
}

void Sample::write_abundances(const std::vector<std::string> &cluster_indicators_to_string, std::string outfile) const {
  // Write relative abundances to a file,
  // outputs to std::cout if outfile is empty.
  const std::vector<double> &abundances = this->group_abundances();

  std::streambuf *buf;
  std::ofstream of;
  if (outfile.empty()) {
    buf = std::cout.rdbuf();
  } else {
    outfile += "_abundances.txt";
    of.open(outfile);
    buf = of.rdbuf();
  }
  std::ostream out(buf);
  out << "#mSWEEP_version:" << '\t' << _BUILD_VERSION << '\n';
  out << "#total_hits:" << '\t' << this->counts_total << '\n';
  out << "#c_id" << '\t' << "mean_theta" << '\n';
  for (size_t i = 0; i < abundances.size(); ++i) {
    out << cluster_indicators_to_string[i] << '\t' << abundances[i] << '\n';
  }
  if (!outfile.empty()) {
    of.close();
  }
}

void Sample::write_bootstrap(const std::vector<std::string> &cluster_indicators_to_string, std::string outfile, unsigned iters) {
  // Write relative abundances to a file,
  // outputs to std::cout if outfile is empty.
  std::streambuf *buf;
  std::ofstream of;
  if (outfile.empty()) {
    buf = std::cout.rdbuf();
  } else {
    outfile += "_abundances.txt";
    of.open(outfile);
    buf = of.rdbuf();
  }
  std::ostream out(buf);
  out << "#c_id" << '\t' << "mean_theta" << '\t' << "abundances" << '\t' << "bootstrap_abundances" << '\n';
  out << "#total_hits:" << '\t' << this->counts_total << '\n';

  for (size_t i = 0; i < cluster_indicators_to_string.size(); ++i) {
    out << cluster_indicators_to_string[i] << '\t';
    for (unsigned j = 0; j < iters; ++j) {
      out << this->bootstrap_abundances.at(j).at(i) << (j == iters - 1 ? '\n' : '\t');
    }
  }
  out << std::endl;
  if (!outfile.empty()) {
    of.close();
  }
}

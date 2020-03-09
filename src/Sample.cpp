#include "Sample.hpp"

#include "likelihood.hpp"
#include "version.h"

void Sample::process_aln() {
  cell_id = "";
  m_num_ecs = pseudos.size();
  m_num_refs = pseudos.n_targets();
  log_ec_counts.resize(m_num_ecs, 0.0);
  uint32_t aln_counts_total = 0;
#pragma omp parallel for schedule(static) reduction(+:aln_counts_total)
  for (uint32_t i = 0; i < m_num_ecs; ++i) {
    log_ec_counts[i] = std::log(pseudos.ec_counts[i]);
    aln_counts_total += pseudos.ec_counts[i];
  }
  counts_total = aln_counts_total;
  //  pseudos.ec_counts.clear();
}

void Sample::read_themisto(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &strands) {
  ReadThemisto(mode, n_refs, strands, &pseudos);
  process_aln();
}

void Sample::read_kallisto(const uint32_t n_refs, std::istream &ec_file, std::istream &tsv_file) {
  ReadKallisto(n_refs, ec_file, tsv_file, &pseudos);
  process_aln();
}

std::vector<double> Sample::group_abundances() const {
  // Calculate the relative abundances of the
  // reference groups from the ec_probs matrix
  std::vector<double> thetas(this->ec_probs.get_rows(), 0.0);
  for (uint32_t i = 0; i < this->ec_probs.get_rows(); ++i) {
    for (uint32_t j = 0; j < this->ec_probs.get_cols(); ++j) {
      thetas[i] += std::exp(this->ec_probs(i, j) + this->log_ec_counts[j]);
    }
    thetas[i] /= this->counts_total;
  }
  return thetas;
}

std::vector<uint16_t> Sample::group_counts(const std::vector<uint32_t> indicators, const uint32_t ec_id, const uint32_t n_groups) const {
  std::vector<uint16_t> read_hitcounts(n_groups);
  for (uint32_t j = 0; j < m_num_refs; ++j) {
    read_hitcounts[indicators[j]] += pseudos.ec_configs[ec_id][j];
  }
  return read_hitcounts;
}

void Sample::write_probabilities(const std::vector<std::string> &cluster_indicators_to_string, const bool gzip_probs, std::ostream &of) const {
  // Write the probability matrix to a file.
  if (of.good()) {
    of << "ec_id" << ',';
    for (uint32_t i = 0; i < this->ec_probs.get_rows(); ++i) {
      of << cluster_indicators_to_string[i];
      of << (i < this->ec_probs.get_rows() - 1 ? ',' : '\n');
    }
    for (uint32_t i = 0; i < this->ec_probs.get_cols(); ++i) {
      of << pseudos.ec_ids[i] << ',';
      for (uint32_t j = 0; j < this->ec_probs.get_rows(); ++j) {
	of << std::exp(this->ec_probs(j, i));
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

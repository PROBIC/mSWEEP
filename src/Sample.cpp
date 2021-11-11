#include "Sample.hpp"

#include <cmath>
#include <sstream>

#include "rcgpar.hpp"

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
}

std::vector<uint16_t> Sample::group_counts(const std::vector<uint32_t> indicators, const uint32_t ec_id, const uint32_t n_groups) const {
  std::vector<uint16_t> read_hitcounts(n_groups);
  for (uint32_t j = 0; j < m_num_refs; ++j) {
    read_hitcounts[indicators[j]] += pseudos.ec_configs[ec_id][j];
  }
  return read_hitcounts;
}

void Sample::write_probabilities(const std::vector<std::string> &cluster_indicators_to_string, std::ostream &of) const {
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
  const std::vector<double> &abundances = rcgpar::mixture_components(this->ec_probs, this->log_ec_counts);

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
  out << "#mSWEEP_version:" << '\t' << MSWEEP_BUILD_VERSION << '\n';
  out << "#total_hits:" << '\t' << this->counts_total << '\n';
  out << "#c_id" << '\t' << "mean_theta" << '\n';
  for (size_t i = 0; i < abundances.size(); ++i) {
    out << cluster_indicators_to_string[i] << '\t' << abundances[i] << '\n';
  }
  if (!outfile.empty()) {
    of.close();
  }
}

void Sample::write_likelihood(const bool gzip_output, const uint32_t n_groups, std::string outfile) const {
  // Write likelihoods to a file

  std::streambuf *buf;
  std::unique_ptr<std::ostream> of;
  if (outfile.empty()) {
    buf = std::cout.rdbuf();
  } else {
    outfile += "_likelihoods.txt";
    switch(gzip_output) { // might want to use other compressions in the future
    case 1:
      outfile += ".gz";
      of = std::unique_ptr<std::ostream>(new bxz::ofstream(outfile));
      break;
    default:
      of = std::unique_ptr<std::ostream>(new std::ofstream(outfile));
      break;
    }
    buf = of->rdbuf();
  }
  std::ostream out(buf);

  for (uint32_t i = 0; i < this->m_num_ecs; ++i){
    uint32_t ec_hit_count = std::round(std::exp(this->log_ec_counts[i]));
    out << ec_hit_count << '\t';
      for (uint32_t j = 0; j < n_groups; ++j) {
	out << this->ll_mat(j, i);
	out << (j == n_groups - 1 ? '\n' : '\t');
      }
  }
  if (!outfile.empty()) {
    of->flush();
  }
}

void Sample::write_likelihood_bitseq(const bool gzip_output, const uint32_t n_groups, std::string outfile) const {
  // Write likelihoods to a file
  // *Note*: will write in BitSeq format!
  // Use Sample::write_likelihoods if tab-separated matrix format is needed.

  std::streambuf *buf;
  std::unique_ptr<std::ostream> of;
  if (outfile.empty()) {
    buf = std::cout.rdbuf();
  } else {
    outfile += "_bitseq_likelihoods.txt";
    switch(gzip_output) { // might want to use other compressions in the future
    case 1:
      outfile += ".gz";
      of = std::unique_ptr<std::ostream>(new bxz::ofstream(outfile));
      break;
    default:
      of = std::unique_ptr<std::ostream>(new std::ofstream(outfile));
      break;
    }
    buf = of->rdbuf();
  }
  std::ostream out(buf);
  out << "# Ntotal " << this->counts_total << '\n';
  out << "# Nmap " << this->counts_total << '\n';
  out << "# M " << n_groups << '\n';
  out << "# LOGFORMAT (probabilities saved on log scale.)" << '\n';
  out << "# r_name num_alignments (tr_id prob )^*{num_alignments}" << '\n';

  uint32_t read_id = 1;
  for (uint32_t i = 0; i < this->m_num_ecs; ++i) {
    uint32_t ec_hit_count = std::round(std::exp(this->log_ec_counts[i]));
    for (uint32_t k = 0; k < ec_hit_count; ++k) {
      out << read_id << ' ';
      out << n_groups + 1 << ' ';
      for (uint32_t j = 0; j < n_groups; ++j) {
	out << j + 1 << ' ' << this->ll_mat(j, i) << ' ';
      }
      out << 0 << ' ' << "-10000.00" << '\n';
      ++read_id;
    }
  }
  if (!outfile.empty()) {
    of->flush();
  }
}

void Sample::CalcLikelihood(const Grouping &grouping, const double bb_constants[2], const std::vector<uint32_t> &group_indicators, const bool cleanup) {
  uint32_t n_groups = grouping.get_n_groups();
  likelihood_array_mat(grouping, group_indicators, bb_constants, *this);
  if (cleanup) {
    // If estimating with only 1 grouping free the memory used by the configs
    clear_configs();
  }
}

void Sample::ReadLikelihood(const Grouping &grouping, std::istream &infile) {
  uint32_t n_groups = grouping.get_n_groups();

  std::vector<std::vector<double>> likelihoods(n_groups, std::vector<double>());
  this->m_num_ecs = 0;
  this->counts_total = 0;
  this->pseudos = KallistoAlignment();

  if (infile.good()) {
    std::string newline;
    uint32_t line_nr = 0;
    while (std::getline(infile, newline)) {
      this->pseudos.ec_ids.emplace_back(line_nr);
      ++line_nr;
      std::string part;
      std::stringstream partition(newline);
      bool ec_count_col = true;
      uint32_t group_id = 0;
      while (std::getline(partition, part, '\t')) {
	if (ec_count_col) {
	  uint32_t ec_count = std::stol(part);
	  this->pseudos.ec_counts.emplace_back(ec_count);
	  this->counts_total += ec_count;
	  this->log_ec_counts.emplace_back(std::log(ec_count));
	  ec_count_col = false;
	} else {
	  likelihoods[group_id].emplace_back(std::stod(part));
	  ++group_id;
	}
      }
    }
    this->m_num_ecs = line_nr;
  } else {
    throw std::runtime_error("Could not read from the likelihoods file.");
  }
  this->ll_mat = rcgpar::Matrix<double>(likelihoods);
}

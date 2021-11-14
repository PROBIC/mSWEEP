#include "Sample.hpp"

#include <cmath>
#include <sstream>

#include "rcgpar.hpp"

#include "version.h"

void Sample::process_aln(const bool bootstrap_mode) {
  cell_id = "";
  m_num_ecs = pseudos.ec_counts.size();
  log_ec_counts.resize(m_num_ecs, 0.0);
  uint32_t aln_counts_total = 0;
#pragma omp parallel for schedule(static) reduction(+:aln_counts_total)
  for (uint32_t i = 0; i < m_num_ecs; ++i) {
    log_ec_counts[i] = std::log(pseudos.ec_counts[i]);
    aln_counts_total += pseudos.ec_counts[i];
  }
  counts_total = aln_counts_total;

  if (!bootstrap_mode) {
    // EC counts aren't needed when not bootstrapping.
    this->pseudos.ec_counts.clear();
    this->pseudos.ec_counts.shrink_to_fit();
  }
}

std::vector<uint16_t> Sample::group_counts(const std::vector<uint32_t> indicators,
					   const uint32_t ec_id, const uint32_t n_groups) const {
  std::vector<uint16_t> read_hitcounts(n_groups);
  uint32_t m_num_refs = this->pseudos.n_targets();
  for (uint32_t j = 0; j < m_num_refs; ++j) {
    read_hitcounts[indicators[j]] += pseudos.ec_configs[ec_id][j];
  }
  return read_hitcounts;
}

void Sample::write_probabilities(const std::vector<std::string> &cluster_indicators_to_string,
				 std::ostream &of) const {
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
    of << std::endl;
    of.flush();
  } else {
    throw std::runtime_error("Can't write to probs file.");
  }
}

void Sample::write_abundances(const std::vector<std::string> &cluster_indicators_to_string, std::ostream &of) const {
  // Write relative abundances to &of,
  if (of.good()) {
    of << "#mSWEEP_version:" << '\t' << MSWEEP_BUILD_VERSION << '\n';
    of << "#total_hits:" << '\t' << this->counts_total << '\n';
    of << "#c_id" << '\t' << "mean_theta" << '\n';
    for (size_t i = 0; i < relative_abundances.size(); ++i) {
      of << cluster_indicators_to_string[i] << '\t' << relative_abundances[i] << '\n';
    }
    of.flush();
  } else {
    throw std::runtime_error("Can't write to abundances file.");
  }
}

void Sample::write_likelihood(const bool gzip_output, const uint32_t n_groups, std::string outfile) const {
  // Write likelihoods to a file

  std::streambuf *buf;
  cxxio::Out of;
  if (outfile.empty()) {
    buf = std::cout.rdbuf();
  } else {
    outfile += "_likelihoods.txt";
    if (gzip_output) {
      outfile += ".gz";
      of.open_compressed(outfile);
    } else {
      of.open(outfile);
    }
    buf = of.stream().rdbuf();
  }
  std::ostream out(buf);

  if (out.good()) {
    for (uint32_t i = 0; i < this->m_num_ecs; ++i){
      uint32_t ec_hit_count = std::round(std::exp(this->log_ec_counts[i]));
      out << ec_hit_count << '\t';
      for (uint32_t j = 0; j < n_groups; ++j) {
	out << this->ll_mat(j, i);
	out << (j == n_groups - 1 ? '\n' : '\t');
      }
    }
    out.flush();
  } else {
    throw std::runtime_error("Can't write to likelihoods file.");
  }
}

void Sample::write_likelihood_bitseq(const bool gzip_output, const uint32_t n_groups, std::string outfile) const {
  // Write likelihoods to a file
  // *Note*: will write in BitSeq format!
  // Use Sample::write_likelihoods if tab-separated matrix format is needed.

  std::streambuf *buf;
  cxxio::Out of;
  if (outfile.empty()) {
    buf = std::cout.rdbuf();
  } else {
    outfile += "_bitseq_likelihoods.txt";
    if (gzip_output) {
      outfile += ".gz";
      of.open_compressed(outfile);
    } else {
      of.open(outfile);
    }
    buf = of.stream().rdbuf();
  }
  std::ostream out(buf);
  if (out.good()) {
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
    out.flush();
  } else {
    throw std::runtime_error("Can't write to likelihoods file (bitseq format).");
  }
}

void Sample::read_likelihood(const Grouping &grouping, std::istream &infile) {
  uint32_t n_groups = grouping.get_n_groups();

  std::vector<std::vector<double>> likelihoods(n_groups, std::vector<double>());
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
	  ec_count_col = false;
	} else {
	  likelihoods[group_id].emplace_back(std::stod(part));
	  ++group_id;
	}
      }
    }
  } else {
    throw std::runtime_error("Could not read from the likelihoods file.");
  }
  this->ll_mat = rcgpar::Matrix<double>(likelihoods);
}

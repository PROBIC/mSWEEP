#include "Sample.hpp"

#include <cmath>
#include <sstream>

#include "rcgpar.hpp"

#include "version.h"

void Sample::process_aln(const telescope::GroupedAlignment &pseudos, const bool bootstrap_mode) {
  cell_id = "";
  m_num_ecs = pseudos.n_ecs();
  log_ec_counts.resize(m_num_ecs, 0.0);
  uint32_t aln_counts_total = 0;
#pragma omp parallel for schedule(static) reduction(+:aln_counts_total)
  for (uint32_t i = 0; i < m_num_ecs; ++i) {
    log_ec_counts[i] = std::log(pseudos.reads_in_ec(i));
    aln_counts_total += pseudos.reads_in_ec(i);
  }
  counts_total = aln_counts_total;

  if (!bootstrap_mode) {
    // EC counts aren't needed when not bootstrapping.
    //this->pseudos.clear_counts();
  }
}

// std::vector<uint16_t> Sample::group_counts(const std::vector<uint32_t> indicators,
// 					   const uint32_t ec_id, const uint32_t n_groups) const {
//   std::vector<uint16_t> read_hitcounts(n_groups);
//   uint32_t m_num_refs = this->pseudos.n_targets();
//   for (uint32_t j = 0; j < m_num_refs; ++j) {
//     read_hitcounts[indicators[j]] += pseudos(ec_id, j);
//   }
//   return read_hitcounts;
// }

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
      of << i << ',';
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

#include "Sample.hpp"

#include "likelihood.hpp"

void BootstrapSample::init_bootstrap(Grouping &grouping) {
  ec_distribution = std::discrete_distribution<uint32_t>(pseudos.aln.ec_counts.begin(), pseudos.aln.ec_counts.end());
  ll_mat = likelihood_array_mat(*this, grouping);
}

void BootstrapSample::resample_counts(std::mt19937_64 &generator) {
  std::vector<uint32_t> tmp_counts(m_num_ecs);
  for (uint32_t i = 0; i < counts_total; ++i) {
    uint32_t ec_id = ec_distribution(generator);
    tmp_counts[ec_id] += 1;
  }
#pragma omp parallel for schedule(static)
  for (uint32_t i = 0; i < m_num_ecs; ++i) {
    log_ec_counts[i] = std::log(tmp_counts[i]);
  }
}

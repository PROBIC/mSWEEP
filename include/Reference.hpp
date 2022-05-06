#ifndef MSWEEP_REFERENCE_HPP
#define MSWEEP_REFERENCE_HPP

#include <cstddef>
#include <vector>
#include <string>
#include <fstream>

#include "cxxio.hpp"

#include "Grouping.hpp"

class Reference {
private:
  uint32_t n_refs = 0;
  uint16_t n_groupings = 0;

  std::vector<std::vector<uint32_t>> groups_indicators;
  std::vector<Grouping> groupings;

  void add_sequence(const std::string &seq_name, const uint16_t grouping_id);

public:
  void read_from_file(std::istream &indicator_file, const char delimiter = '\t');
  void match_with_fasta(const char delimiter, std::istream &groups_file, std::istream &fasta_file);

  // Check that the reference has same number of sequences as ...
  // ... the themisto index.
  void verify_themisto_index(cxxio::In &themisto_index) const;
  // ... the kallisto pseudoalignment.
  void verify_kallisto_alignment(std::istream &kallisto_run_info) const;

  // Getters
  const Grouping& get_grouping(const uint16_t grouping_id) const { return this->groupings[grouping_id]; };
  const std::vector<uint32_t>& get_group_indicators(const uint16_t grouping_id) const { return this->groups_indicators[grouping_id]; };
  uint32_t get_n_refs() const { return this->n_refs; };
  uint16_t get_n_groupings() const { return this->n_groupings; };
};

#endif

#ifndef MSWEEP_REFERENCE_HPP
#define MSWEEP_REFERENCE_HPP

#include <vector>
#include <fstream>
#include <string>

#include "file.hpp"

#include "Grouping.hpp"

class Reference {
public:
  uint32_t n_refs = 0;
  std::vector<std::vector<uint32_t>> groups_indicators;
  std::vector<Grouping> groupings;

  void read_from_file(std::istream &indicator_file, const char delimiter = '\t');
  void match_with_fasta(const char delimiter, std::istream &groups_file, std::istream &fasta_file);
  void add_sequence(const std::string &seq_name, const uint8_t grouping_id);

  void verify_themisto_index(File::In &themisto_index) const;
  void verify_kallisto_alignment(std::istream &kallisto_run_info) const;
};

#endif

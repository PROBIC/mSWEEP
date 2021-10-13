#ifndef MSWEEP_REFERENCE_HPP
#define MSWEEP_REFERENCE_HPP

#include <fstream>

#include "file.hpp"

#include "Grouping.hpp"

class Reference {
public:
  uint32_t n_refs = 0;
  std::vector<uint32_t> group_indicators;
  Grouping grouping;

  void add_sequence(const std::string &seq_name);

  void verify(std::istream &infile) const;
  void verify(File::In &infile) const;
private:
  void verify_themisto_index(File::In &themisto_index) const;
  void verify_kallisto_alignment(std::istream &kallisto_run_info) const;
};

void ReadClusterIndicators(std::istream &indicator_file, Reference &reference);
void MatchClusterIndicators(const char delim, std::istream &groups, std::istream &fasta, Reference &reference);

#endif

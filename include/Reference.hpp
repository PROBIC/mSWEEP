#ifndef MSWEEP_REFERENCE_HPP
#define MSWEEP_REFERENCE_HPP

#include <fstream>

#include "Grouping.hpp"

class Reference {
public:
  Grouping grouping;
  uint32_t n_refs;

  void verify(const bool themisto_mode, std::istream &infile) const;
private:
  void verify_themisto_index(std::istream &themisto_index) const;
  void verify_kallisto_alignment(std::istream &kallisto_run_info) const;
};

void ReadClusterIndicators(std::istream &indicator_file, Reference &reference);
void MatchClusterIndicators(const char delim, std::istream &groups, std::istream &fasta, Reference &reference);

#endif

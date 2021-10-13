#ifndef MSWEEP_REFERENCE_HPP
#define MSWEEP_REFERENCE_HPP

#include <vector>
#include <array>
#include <string>
#include <fstream>

#include "Grouping.hpp"

class Reference {
public:
  Grouping grouping;
  std::vector<std::string> group_names;
  uint32_t n_refs;

  void verify() const;
};

void ReadClusterIndicators(std::istream &indicator_file, Reference &reference);
void MatchClusterIndicators(const char delim, std::istream &groups, std::istream &fasta, Reference &reference);

#endif

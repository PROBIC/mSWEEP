#include "Reference.hpp"

#include <unordered_map>
#include <exception>

#include "tools/matchfasta.hpp"

void Reference::verify() const {
  if (this->n_refs == 0) {
    throw std::runtime_error("The grouping contains 0 reference sequences");
  }
}

void AddToReference(const std::string &indicator_s, std::unordered_map<std::string, unsigned> &str_to_int, Reference &reference) {
  if (str_to_int.find(indicator_s) == str_to_int.end()) {
    str_to_int[indicator_s] = str_to_int.size();
    reference.group_names.emplace_back(indicator_s);
    reference.grouping.sizes.emplace_back(0);
  }
  ++reference.grouping.sizes[str_to_int[indicator_s]];
  reference.grouping.indicators.emplace_back(str_to_int[indicator_s]);
}

void ReadClusterIndicators(std::istream &indicator_file, Reference &reference) {
  std::unordered_map<std::string, unsigned> str_to_int;

  if (indicator_file.good()) {
    std::string indicator_s;
    while (getline(indicator_file, indicator_s)) {
      AddToReference(indicator_s, str_to_int, reference);
    }
  } else {
    throw std::runtime_error("Could not read cluster indicators.");
  }

  reference.n_refs = reference.grouping.indicators.size();
  reference.grouping.n_groups = str_to_int.size();
}

void MatchClusterIndicators(const char delim, std::istream &groups, std::istream &fasta, Reference &reference) {
  std::unordered_map<std::string, unsigned> str_to_int;
  std::vector<std::string> groups_in_fasta;
  try {
    mSWEEP::tools::matchfasta(groups, fasta, delim, &groups_in_fasta);
  } catch (std::exception &e) {
    throw std::runtime_error("Matching the group indicators to the fasta file failed, is the --groups-list delimiter correct?");
  }

  for (uint32_t i = 0; i < groups_in_fasta.size(); ++i) {
    AddToReference(groups_in_fasta[i], str_to_int, reference);
  }

  reference.n_refs = reference.grouping.indicators.size();
  reference.grouping.n_groups = str_to_int.size();
}

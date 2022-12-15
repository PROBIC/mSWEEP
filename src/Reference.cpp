#include "Reference.hpp"

#include <unordered_map>
#include <exception>
#include <sstream>

#include "tools/matchfasta.hpp"

void Reference::add_sequence(const std::string &indicator_s, const uint16_t grouping_id) {
  this->groupings[grouping_id].add_sequence(indicator_s);
  if (grouping_id == 0) {
    this->n_refs += 1;
  }
  this->groups_indicators[grouping_id].emplace_back(this->groupings[grouping_id].get_id(indicator_s));
}

void Reference::read_from_file(std::istream &indicator_file, const char delimiter) {
  if (indicator_file.good()) {
    std::string indicator_s;
    while (std::getline(indicator_file, indicator_s)) {
      std::stringstream indicators(indicator_s);
      std::string indicator;
      uint16_t grouping_id = 0;
      while (std::getline(indicators, indicator, delimiter)) {
	if (grouping_id >= this->n_groupings) {
	  this->groupings.emplace_back(Grouping());
	  this->groups_indicators.emplace_back(std::vector<uint32_t>());
	  ++this->n_groupings;
	}
	this->add_sequence(indicator, grouping_id);
	++grouping_id;
      }
    }
  } else {
    throw std::runtime_error("Could not read cluster indicators.");
  }
  if (this->n_refs == 0) {
    throw std::runtime_error("The grouping contains 0 reference sequences");
  }
}

void Reference::match_with_fasta(const char delim, std::istream &groups_file, std::istream &fasta_file) {
  std::vector<std::vector<std::string>> groups_in_fasta;
  try {
    mSWEEP::tools::matchfasta(groups_file, fasta_file, delim, &groups_in_fasta);
  } catch (std::exception &e) {
    throw std::runtime_error("Matching the group indicators to the fasta file failed, is the --groups-delimiter argument correct?");
  }

  this->n_groupings = groups_in_fasta[0].size();
  this->groupings = std::vector<Grouping>(this->n_groupings);
  this->groups_indicators = std::vector<std::vector<uint32_t>>(this->n_groupings);

  for (uint32_t i = 0; i < groups_in_fasta.size(); ++i) {
    for (uint16_t j = 0; j < this->n_groupings; ++j) {
      this->add_sequence(groups_in_fasta[i][j], j);
    }
  }
}

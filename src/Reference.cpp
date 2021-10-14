#include "Reference.hpp"

#include <unordered_map>
#include <exception>
#include <sstream>

#include "tools/matchfasta.hpp"

void Reference::verify_themisto_index(File::In &themisto_index) const {
  uint32_t lines_in_grouping = themisto_index.count_lines<uint32_t>();
  if (lines_in_grouping > this->n_refs) {
    throw std::runtime_error("pseudoalignment has more reference sequences than the grouping.");
  } else if (lines_in_grouping < this->n_refs) {
    throw std::runtime_error("grouping has more reference sequences than the pseudoalignment.");
  }
}

void Reference::verify_kallisto_alignment(std::istream &run_info) const {
  // Get the number of reference sequences in the pseudoalignment
  // contained in the 'n_targets' variable in run_info.json file.
  short unsigned line_nr = 0; // number of reference seqs is on line 2 (kallisto v0.43)
  if (run_info.good()) {
    std::string line;
    while (getline(run_info, line)) {
      if (line_nr == 0) {
	++line_nr;
      } else {
	std::string part;
	std::stringstream partition(line);
	unsigned n_targets = 0;
	while (getline(partition, part, ':')) {
	  if (n_targets == 0) {
	    ++n_targets;
	  } else {
	    part.pop_back(); // the number ends in a ','; get rid of it.
	    unsigned n_targets = std::stoi(part);
	    if (n_targets > this->n_refs) {
	      throw std::runtime_error("pseudoalignment has more reference sequences than the grouping.");
	    } else if (n_targets < this->n_refs) {
	      throw std::runtime_error("grouping has more reference sequences than the pseudoalignment.");
	    }
	    return;
	  }
	}
      }
    }
  } else {
    throw std::runtime_error("Could not read run_info.json found.");
  }
}

void Reference::add_sequence(const std::string &indicator_s) {
  this->grouping.add_group(indicator_s);
  this->grouping.add_sequence(indicator_s);
  this->n_refs += 1;
  this->group_indicators.emplace_back(this->grouping.name_to_id[indicator_s]);
}

void Reference::read_from_file(std::istream &indicator_file) {
  if (indicator_file.good()) {
    std::string indicator_s;
    while (getline(indicator_file, indicator_s)) {
      this->add_sequence(indicator_s);
    }
  } else {
    throw std::runtime_error("Could not read cluster indicators.");
  }
  if (this->n_refs == 0) {
    throw std::runtime_error("The grouping contains 0 reference sequences");
  }
}

void Reference::match_with_fasta(const char delim, std::istream &groups_file, std::istream &fasta_file) {
  std::vector<std::string> groups_in_fasta;
  try {
    mSWEEP::tools::matchfasta(groups_file, fasta_file, delim, &groups_in_fasta);
  } catch (std::exception &e) {
    throw std::runtime_error("Matching the group indicators to the fasta file failed, is the --groups-list delimiter correct?");
  }

  for (uint32_t i = 0; i < groups_in_fasta.size(); ++i) {
    this->add_sequence(groups_in_fasta[i]);
  }
}

#include "matchfasta.hpp"

#include <sstream>
#include <utility>
#include <algorithm>
#include <unordered_map>

namespace mSWEEP {
namespace tools {

void read_groups(const char delim, std::istream &groups, std::unordered_map<std::string, std::vector<std::string>> *seq_to_group) {
  std::string line;
  while(getline(groups, line, '\n')) {
    std::string seqname;
    std::string group;
    std::stringstream parts(line);
    getline(parts, seqname, delim);
    if (seq_to_group->find(seqname) == seq_to_group->end()) {
      std::getline(parts, group, delim);
      std::vector<std::string> groups(1, group);
      seq_to_group->insert(std::make_pair(seqname, groups));
      while (std::getline(parts, group, delim)) {
	(*seq_to_group)[seqname].emplace_back(group);
      }
    }
  }
}
void read_fasta(const std::unordered_map<std::string, std::vector<std::string>> &seq_to_group, std::istream &fasta, std::vector<std::vector<std::string>> *groups_in_fasta) {
  std::string line;
  size_t n_groupings = seq_to_group.begin()->second.size();
  while(getline(fasta, line, '\n')) {
    if (line.at(0) == '>') { // todo check if seqence exists in the unordered_map
      const std::string seq_name = line.erase(0, 1);
      groups_in_fasta->emplace_back(std::vector<std::string>());
      for (size_t i = 0; i < n_groupings; ++i) {
	groups_in_fasta->back().emplace_back(seq_to_group.at(seq_name)[i]);
      }
    }
  }
}
void matchfasta(std::istream &groups, std::istream &fasta, const char delim, std::vector<std::vector<std::string>> *groups_in_fasta) {
  std::unordered_map<std::string, std::vector<std::string>> seq_to_group;
  read_groups(delim, groups, &seq_to_group);
  read_fasta(seq_to_group, fasta, groups_in_fasta);
}
}
}

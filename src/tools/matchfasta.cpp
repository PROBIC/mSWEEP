#include "matchfasta.hpp"

#include <sstream>
#include <utility>
#include <algorithm>
#include <unordered_map>

namespace mSWEEP {
namespace tools {

void read_groups(std::istream &groups, const char delim, std::unordered_map<std::string, std::string> *seq_to_group) {
  std::string line;
  while(getline(groups, line, '\n')) {
    std::string seqname;
    std::string group;
    std::stringstream parts(line);
    getline(parts, seqname, delim);
    if (seq_to_group->find(seqname) == seq_to_group->end()) {
      getline(parts, group, delim);
      seq_to_group->insert(std::make_pair(seqname, group));
    }
  }
}
void read_fasta(std::istream &fasta, const std::unordered_map<std::string, std::string> &seq_to_group, std::vector<std::string> *groups_in_fasta) {
  std::string line;
  while(getline(fasta, line, '\n')) {
    if (line.at(0) == '>') { // todo check if seqence exists in the unordered_map
      groups_in_fasta->push_back(seq_to_group.at(line.erase(0, 1)));
    }
  }
}
void matchfasta(std::istream &groups, std::istream &fasta, const char delim, std::vector<std::string> *groups_in_fasta) {
  std::unordered_map<std::string, std::string> seq_to_group;
  read_groups(groups, delim, &seq_to_group);
  read_fasta(fasta, seq_to_group, groups_in_fasta);
}
}
}

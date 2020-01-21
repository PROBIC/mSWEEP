#include <string>
#include <exception>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <utility>
#include <vector>
#include <algorithm>

namespace tools {
namespace exceptions {
struct tools : public std::exception {
  std::string msg;
  tools(tools&&) = default;
  tools(std::string m) : msg(m) {}
  const char* what() const noexcept override { return msg.c_str(); }
};
struct ifstream : public tools {
  ifstream(std::string m) : tools("Cannot read from file " + m) {}
};
struct notfound : public tools {
  notfound(std::string m) : tools("Argument " + m + " not given.") {};
};
}

class Args {
 private:
  char* GetCmdOption(char **begin, char **end, const std::string &option) const {
    char **it = std::find(begin, end, option);
    return ((it != end && ++it != end) ? *it : 0);
  }
  bool CmdOptionPresent(char **begin, char **end, const std::string &option) const {
    return (std::find(begin, end, option) != end);
  }
  void print_help() {
    std::cerr << "Usage: matchfasta --fasta <fasta file> --groups <indictators table> -d <delimiter (default:tab)>\n"
      << "Match first column of --groups to the lines in --fasta, and print second column of --groups in that order.\n\n"
      << "Options:\n"
      << "\t--fasta <fastaFile>\n"
      << "\tFasta file with sequence names defined by lines starting with '>'\n"
      << "\t--groups <groupsFile>\n"
      << "\tTable containing sequence names in the first column and group indicators in the second.\n"
      << "\t-d <delimiterChar>>\n"
      << "\tDelimiter character for --groups (default: 1 tab)\n" << std::endl;
  }

 public:
  std::ifstream fasta;
  std::ifstream groups;
  char delim = '\t';
  Args() {}
  void parse_args(int argc, char *argv[]) {
    if (CmdOptionPresent(argv, argv+argc, "--help")) {
      print_help();
    }
    std::string fasta_file;
    std::string groups_file;

    if (CmdOptionPresent(argv, argv+argc, "--fasta")) {
      fasta_file = std::string(GetCmdOption(argv, argv+argc, "--fasta"));
      fasta.open(fasta_file);
    } else {
      throw exceptions::notfound("--fasta");
    }
    if (fasta.fail()) {
      throw exceptions::ifstream(fasta_file);
    }

    if (CmdOptionPresent(argv, argv+argc, "--groups")) {
      groups_file = std::string(GetCmdOption(argv, argv+argc, "--groups"));
      groups.open(groups_file);
    } else {
      throw exceptions::notfound("--groups");
    }
    if (groups.fail()) {
      throw exceptions::ifstream(groups_file);
    }

    if (CmdOptionPresent(argv, argv+argc, "-d")) {
      delim = *GetCmdOption(argv, argv+argc, "-d");
    }
  }
};

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
void match_groups(std::istream &groups, std::istream &fasta, const char delim, std::vector<std::string> *groups_in_fasta) {
  std::unordered_map<std::string, std::string> seq_to_group;
  read_groups(groups, delim, &seq_to_group);
  read_fasta(fasta, seq_to_group, groups_in_fasta);
}
}

int main (int argc, char *argv[]) {
  tools::Args args;
  try {
    args.parse_args(argc, argv);
  } catch (std::exception &e) {
    std::cerr << "Parsing arguments failed:\n  "
	      << e.what() << '\n'
	      << "exiting" << std::endl;
    return 1;
  }

  std::vector<std::string> groups_in_fasta;
  try {
    tools::match_groups(args.groups, args.fasta, args.delim, &groups_in_fasta);
  } catch (std::exception &e) {
    std::cerr << "Matching indicators failed :\n  "
	      << e.what() << '\n'
	      << "exiting" << std::endl;
    return 1;
  }

  for (auto val : groups_in_fasta) {
    std::cout << val << '\n';
  }
  std::cout << std::endl;
  return 0;
}

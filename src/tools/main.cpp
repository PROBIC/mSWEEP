#include <exception>
#include <vector>
#include <string>
#include <iostream>

#include "matchfasta.hpp"

int main (int argc, char *argv[]) {
  mSWEEP::tools::Args args;
  try {
    args.parse_args(argc, argv);
  } catch (std::exception &e) {
    std::cerr << "Parsing arguments failed:\n  "
	      << e.what() << '\n'
	      << "exiting" << std::endl;
    return 1;
  }

  std::vector<std::vector<std::string>> groups_in_fasta;
  try {
    mSWEEP::tools::matchfasta(args.groups, args.fasta, args.delim, &groups_in_fasta);
  } catch (std::exception &e) {
    std::cerr << "Matching indicators failed :\n  "
	      << e.what() << '\n'
	      << "exiting" << std::endl;
    return 1;
  }

  size_t n_groups_for_seq = groups_in_fasta.size();
  for (size_t i = 0; i < n_groups_for_seq; ++i) {
    for (auto val : groups_in_fasta[i]) {
      std::cout << val;
      if (i < n_groups_for_seq - 1) {
	std::cout << args.delim;
      }
    }
    std::cout << '\n';
  }
  std::cout << std::endl;
  return 0;
}

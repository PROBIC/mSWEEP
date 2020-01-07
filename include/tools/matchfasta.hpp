#ifndef MSWEEP_TOOLS_MATCHFASTA_HPP
#define MSWEEP_TOOLS_MATCHFASTA_HPP

void match_groups(std::istream &groups, std::istream &fasta, const char delim, std::vector<std::string> *groups_in_fasta);

#endif

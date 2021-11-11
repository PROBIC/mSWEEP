#include "read_pseudoalignment.hpp"

#include <sstream>
#include <exception>
#include <fstream>

#include "bxzstr.hpp"
#include "cxxio.hpp"

std::vector<std::string> ReadCellNames(std::istream &cells_file) {
  Reference reference;
  reference.read_from_file(cells_file);
  return reference.get_grouping(0).get_names();
}

void ReadPseudoalignment(KallistoFiles &kallisto_files, unsigned n_refs, std::vector<std::unique_ptr<Sample>> &batch, bool bootstrap_mode) {
  if (bootstrap_mode) {
    batch.emplace_back(new BootstrapSample());
  } else {
    batch.emplace_back(new Sample());
  }
  ReadKallisto(n_refs, *kallisto_files.ec, *kallisto_files.tsv, &batch.back()->pseudos);
  batch.back()->process_aln();
}

void ReadPseudoalignment(const std::string &tinfile1, const std::string &tinfile2, const std::string &themisto_mode, const bool bootstrap_mode, const unsigned n_refs, std::vector<std::unique_ptr<Sample>> &batch) {
  std::vector<std::istream*> strands(2);
  cxxio::In check_strand_1(tinfile1);
  cxxio::In check_strand_2(tinfile2);
  strands.at(0) = new bxz::ifstream(tinfile1);
  strands.at(1) = new bxz::ifstream(tinfile2);

  if (bootstrap_mode) {
    batch.emplace_back(new BootstrapSample());
  } else {
    batch.emplace_back(new Sample());
  }

  ReadThemisto(get_mode(themisto_mode), n_refs, strands, &batch.back()->pseudos);
  batch.back()->process_aln();

  if (!bootstrap_mode) {
    batch.back()->pseudos.ec_counts.clear();
    batch.back()->pseudos.ec_counts.shrink_to_fit();
  }
}

void ReadLikelihood(const Grouping &grouping, const bool bootstrap_mode, std::istream &infile, std::vector<std::unique_ptr<Sample>> &samples) {
  if (bootstrap_mode) {
    samples.emplace_back(new BootstrapSample());
  } else {
    samples.emplace_back(new Sample());
  }
  samples.back()->ReadLikelihood(grouping, infile);
}

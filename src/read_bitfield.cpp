#include "read_bitfield.hpp"

#include <sstream>
#include <exception>
#include <fstream>

#include "bxzstr.hpp"
#include "file.hpp"

std::vector<std::string> ReadCellNames(std::istream &cells_file) {
  Reference reference;
  ReadClusterIndicators(cells_file, reference);
  return reference.grouping.names;
}

void ReadBitfield(KallistoFiles &kallisto_files, unsigned n_refs, std::vector<std::unique_ptr<Sample>> &batch, Reference &reference, bool bootstrap_mode) {
  if (bootstrap_mode) {
    batch.emplace_back(new BootstrapSample());
  } else {
    batch.emplace_back(new Sample());
  }
  batch.back()->read_kallisto(n_refs, *kallisto_files.ec, *kallisto_files.tsv);
}

void ReadBitfield(const std::string &tinfile1, const std::string &tinfile2, const std::string &themisto_mode, const bool bootstrap_mode, const unsigned n_refs, std::vector<std::unique_ptr<Sample>> &batch) {
  std::vector<std::istream*> strands(2);
  File::In check_strand_1(tinfile1);
  File::In check_strand_2(tinfile2);
  strands.at(0) = new bxz::ifstream(tinfile1);
  strands.at(1) = new bxz::ifstream(tinfile2);

  if (bootstrap_mode) {
    batch.emplace_back(new BootstrapSample());
  } else {
    batch.emplace_back(new Sample());
  }

  batch.back()->read_themisto(get_mode(themisto_mode), n_refs, strands);
}

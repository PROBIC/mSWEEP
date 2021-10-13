#include "read_bitfield.hpp"

#include <sstream>
#include <exception>

#include "bxzstr.hpp"
#include "file.hpp"

void VerifyGrouping(const unsigned n_refs, std::istream &run_info) {
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
	    if (n_targets > n_refs) {
	      throw std::runtime_error("pseudoalignment has more reference sequences than the grouping.");
	    } else if (n_targets < n_refs) {
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

uint32_t CountLines (std::istream &stream) {
  uint32_t n_lines = 0;
  std::string line;
  while (std::getline(stream, line)) {
    n_lines += 1;
  }
  return n_lines;
}

void VerifyThemistoGrouping(const unsigned n_refs, std::istream &themisto_index) {
  uint32_t lines_in_grouping = CountLines(themisto_index);
  if (lines_in_grouping > n_refs) {
    throw std::runtime_error("pseudoalignment has more reference sequences than the grouping.");
  } else if (lines_in_grouping < n_refs) {
    throw std::runtime_error("grouping has more reference sequences than the pseudoalignment.");
  }
}

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

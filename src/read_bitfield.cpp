#include "read_bitfield.hpp"

#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <exception>

void VerifyGrouping(std::istream &run_info, unsigned n_refs) {
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

void ReadClusterIndicators(std::istream &indicator_file, Reference &reference) {
  std::unordered_map<std::string, unsigned> str_to_int;

  if (indicator_file.good()) {
    std::string indicator_s;
    unsigned indicator_i = 0;
    while (getline(indicator_file, indicator_s)) {
      if (str_to_int.find(indicator_s) == str_to_int.end()) {
	str_to_int[indicator_s] = indicator_i;
	reference.group_names.emplace_back(indicator_s);
	reference.grouping.sizes.emplace_back(0);
	++indicator_i;
      }
      ++reference.grouping.sizes[str_to_int[indicator_s]];
      reference.grouping.indicators.emplace_back(str_to_int[indicator_s]);
    }
  } else {
    throw std::runtime_error("Could not read cluster indicators.");
  }

  reference.n_refs = reference.grouping.indicators.size();
  reference.grouping.n_groups = str_to_int.size();
}

std::vector<std::string> ReadCellNames(std::istream &cells_file) {
  Reference reference;
  ReadClusterIndicators(cells_file, reference);
  return reference.group_names;
}

void ReadBitfield(KallistoFiles &kallisto_files, unsigned n_refs, std::vector<Sample> &batch, Reference &reference, bool bootstrap_mode) {
  if (bootstrap_mode) {
    batch.emplace_back(BootstrapSample());
  } else {
    batch.emplace_back(Sample());
  }
  batch.back().read_kallisto(n_refs, *kallisto_files.ec, *kallisto_files.tsv);
}

void ReadBitfield(const std::string &tinfile1, const std::string &tinfile2, const std::string &themisto_mode, const unsigned n_refs, std::vector<Sample> &batch) {
  std::vector<std::istream*> strands(2);
  strands.at(0) = new zstr::ifstream(tinfile1);
  strands.at(1) = new zstr::ifstream(tinfile2);

  batch.emplace_back(Sample());
  batch.back().read_themisto(get_mode(themisto_mode), n_refs, strands);
}

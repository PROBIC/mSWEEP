#include "read_bitfield.hpp"

#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <exception>

#include "telescope/include/telescope.hpp"

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

void ReadBitfield(KallistoFiles &kallisto_files, unsigned n_refs, std::vector<Sample> &batch, Reference &reference) {
  // Reads the alignment file for further analyses.
  // Args:
  //   kallisto_files: kallisto alignment files.
  //   n_refs: number of reference sequences.
  //   batch: alignments for the sample will be appended to this vector
  std::vector<std::string> batch_cells;
  if (kallisto_files.batch_mode) {
    batch_cells = ReadCellNames(*kallisto_files.cells);
  } else {
    batch_cells.emplace_back("sample");
  }

  batch.emplace_back(Sample());
  Sample *current_sample = &batch.back();
  (*current_sample).n_groups = reference.grouping.n_groups;
  (*current_sample).m_num_refs = reference.n_refs;
  if (kallisto_files.tsv->good()) {
    std::string line;
    int cell_id = 0;
    
    while (getline(*kallisto_files.tsv, line)) {
      std::string part;
      std::stringstream partition(line);
      int key = 0;
      int howmany = 0;
      getline(partition, part, '\t');
      key = std::stoi(part);
      getline(partition, part, '\t');
      if (kallisto_files.batch_mode) {
      	int current_cell_id = std::stoi(part);
      	if (current_cell_id != cell_id) {
      	  unsigned num_ecs = (*current_sample).ec_ids.size();
      	  (*current_sample).ec_configs.resize(num_ecs, std::vector<bool>(reference.n_refs, false));
      	  (*current_sample).m_num_ecs = num_ecs;
      	  (*current_sample).cell_id = cell_id;
      	  batch.emplace_back(Sample());
      	  current_sample = &batch.back();
      	  (*current_sample).n_groups = reference.grouping.n_groups;
      	  (*current_sample).m_num_refs = reference.n_refs;
      	  ++cell_id;
      	}
      	getline(partition, part, '\t');
      	howmany = std::stoi(part);
      } else {
	howmany = std::stoi(part);
      }
      if (howmany > 0) {
	(*current_sample).ec_ids.push_back(key);
	(*current_sample).ec_counts.push_back(std::log(howmany));
	(*current_sample).counts_total += howmany;
      }
    }
    unsigned num_ecs = (*current_sample).ec_ids.size();
    (*current_sample).ec_configs.resize(num_ecs, std::vector<bool>(reference.n_refs, false));
    (*current_sample).m_num_ecs = num_ecs;
    (*current_sample).cell_id = cell_id;
  } else {
    throw std::runtime_error(".tsv file not found.");
  }
  if (kallisto_files.ec->good()) {
    std::string line;
    unsigned current_id = 0;
    while (getline(*kallisto_files.ec, line)) {
      std::string part;
      std::stringstream partition(line);
      getline(partition, part, '\t');
      while (getline(partition, part, '\t')) {
	std::string one;
	std::stringstream ones(part);
	while (getline(ones, one, ',')) {
	  unsigned short makeone = std::stoi(one);
	  batch.back().ec_configs[current_id][makeone] = true;
	}
      }
      ++current_id;
    }
  } else {
    throw std::runtime_error(".ec file not found.");
  }
}

void ReadBitfield(const std::string &tinfile1, const std::string &tinfile2, const std::string &themisto_mode, const unsigned n_refs, std::vector<Sample> &batch) {
  std::vector<std::istream*> strands(2);
  strands.at(0) = new zstr::ifstream(tinfile1);
  strands.at(1) = new zstr::ifstream(tinfile2);
  KAlignment kallisto = ReadAlignments(get_mode(themisto_mode), n_refs, &strands);
  batch.emplace_back(Sample(kallisto));
}

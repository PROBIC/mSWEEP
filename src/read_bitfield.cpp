#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <exception>

#include "read_bitfield.hpp"
#include "zstr/zstr.hpp"

void VerifyGrouping(std::string &run_info_file, unsigned n_refs) {
  // Get the number of reference sequences in the pseudoalignment
  // contained in the 'n_targets' variable in run_info.json file.
  short unsigned line_nr = 0; // number of reference seqs is on line 2 (kallisto v0.43)
  std::ifstream run_info(run_info_file);
  if (run_info.is_open()) {
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
    throw std::runtime_error(run_info_file + " not found.");
  }
  run_info.close();
}

void ReadClusterIndicators(std::string &indicator_path, Reference &reference) {
  std::unordered_map<std::string, unsigned> str_to_int;

  zstr::ifstream indicator_file(indicator_path);
  if (indicator_file.good()) {
    std::string indicator_s;
    signed indicator_i = 0;
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
    throw std::runtime_error(indicator_path + " not found.");
  }

  reference.n_refs = reference.grouping.indicators.size();
  reference.grouping.n_groups = str_to_int.size();
}

std::vector<std::string> ReadCellNames(std::string &cells_file) {
  Reference reference;
  ReadClusterIndicators(cells_file, reference);
  return reference.group_names;
}

void ReadBitfield(std::vector<std::string> &kallisto_files, unsigned n_refs, std::vector<Sample> &batch) {
  // Reads the alignment file for further analyses.
  // Args:
  //   kallisto_files: kallisto alignment files.
  //   n_refs: number of reference sequences.
  //   batch: alignments for the sample will be appended to this vector
  bool batch_mode = (kallisto_files.size() == 4); // Batch mode (kallisto v0.44) has four files
  std::vector<std::string> batch_cells;
  if (batch_mode) {
    batch_cells = ReadCellNames(kallisto_files[3]);
  } else {
    batch_cells.emplace_back("sample");
  }
  zstr::ifstream ec_file(kallisto_files[1]);
  zstr::ifstream tsv_file(kallisto_files[2]);

  std::shared_ptr<std::unordered_map<long unsigned, std::vector<bool>>> kallisto_configs;
  kallisto_configs = std::make_shared<std::unordered_map<long unsigned, std::vector<bool>>>();
  std::unordered_set<long unsigned> config_ids;

  if (tsv_file.good()) {
    std::string line;
    unsigned cell_id = 0;
    std::vector<long unsigned> ec_ids;
    std::vector<long unsigned> ec_counts;
    long unsigned counts_total = 0;
    
    while (getline(tsv_file, line)) {
      std::string part;
      std::stringstream partition(line);
      unsigned count = 0;
      long unsigned key = 0;
      long unsigned howmany = 0;
      while (getline(partition, part, '\t')) {
	if (count == 0) {
	  key = std::stoi(part);
	  ++count;
	} else if (count == 1 && batch_mode) {
	  unsigned current_cell_id = std::stoi(part);
	  ++count;
	  if (current_cell_id != cell_id) {
	    batch.emplace_back(Sample(batch_cells[cell_id], ec_ids, ec_counts, counts_total, kallisto_configs));
	    ec_ids.clear();
	    ec_counts.clear();
	    counts_total = 0;
	    ++cell_id;
	  }
	} else {
	  howmany = std::stoi(part);
	}
      }
      if (howmany > 0) {
	config_ids.insert(key);
	ec_ids.push_back(key);
	ec_counts.push_back(howmany);
	counts_total += howmany;
      }
    }
    batch.emplace_back(Sample(batch_cells[cell_id], ec_ids, ec_counts, counts_total, kallisto_configs));
  } else {
    throw std::runtime_error(kallisto_files[2] + " not found.");
  }

  if (ec_file.good()) {
    std::string line;
    while (getline(ec_file, line)) {
      std::string part;
      std::stringstream partition(line);
      bool firstel = true;
      long unsigned key = 0;
      std::vector<bool> config(n_refs, 0);
      bool lookup;
      while (getline(partition, part, '\t')) {
	if (firstel) {
	  key = std::stoi(part);
	  firstel = false;
	  lookup = config_ids.find(key) != config_ids.end();
	} else if (lookup) {
	  std::string one;
	  std::stringstream ones(part);
	  while (getline(ones, one, ',')) {
	    unsigned makeone = std::stoi(one);
	    config[makeone] = 1;
	  }
	}
      }
      if (lookup) {
	(*kallisto_configs)[key] = config;
      }
    }
  } else {
    throw std::runtime_error(kallisto_files[1] + " not found.");
  }
}

#include "mSWEEP.hpp"

#include <iostream>
#include <exception>
#include <fstream>

#include "cxxio.hpp"
#include "rcgpar.hpp"

#include "likelihood.hpp"
#include "version.h"

void ReadGroupIndicators(const std::string &indicators_path, Reference *reference) {
  cxxio::In indicators_file(indicators_path);
  reference->read_from_file(indicators_file.stream(), '\t'); // TODO: take delimiter as argument.
}

telescope::GroupedAlignment ReadPseudoalignments(const std::vector<std::string> &alignment_paths, const std::string &themisto_merge_mode, const Reference &reference) {
  size_t n_files = alignment_paths.size();
  std::vector<cxxio::In> infiles;
  infiles.reserve(n_files);
  std::vector<std::istream*> strands(n_files);
  if (n_files > 0) {
    for (size_t i = 0; i < n_files; ++i) {
      infiles.emplace_back(cxxio::In(alignment_paths[i]));
      strands[i] = &infiles[i].stream();
    }
  } else {
    strands.emplace_back(&std::cin);
  }
  return telescope::read::ThemistoGrouped(telescope::get_mode(themisto_merge_mode), reference.get_n_refs(), reference.get_group_indicators(0), strands);
}

seamat::DenseMatrix<double> ReadLikelihoodFromFile(const std::string &likelihood_path, const Reference &reference, std::ostream &log, std::vector<double> *log_ec_counts) {
  log << "  reading likelihoods from file" << '\n';
  if (reference.get_n_groupings() > 1) {
    throw std::runtime_error("Using more than one grouping with --read-likelihood is not yet implemented.");
  }
  cxxio::In infile(likelihood_path);

  // Have to read the likelihoods into a temporary because num_ecs is not known
  uint32_t n_groups = reference.get_grouping(0).get_n_groups();
  std::vector<std::vector<double>> likelihoods(n_groups, std::vector<double>());

  if (infile.stream().good()) {
    std::string newline;
    uint32_t line_nr = 0;
    while (std::getline(infile.stream(), newline)) {
      ++line_nr;
      std::string part;
      std::stringstream partition(newline);
      bool ec_count_col = true;
      uint32_t group_id = 0;
      while (std::getline(partition, part, '\t')) {
	if (ec_count_col) {
	  uint32_t ec_count = std::stol(part);
	  log_ec_counts->emplace_back(std::log(ec_count));
	  ec_count_col = false;
	} else {
	  likelihoods[group_id].emplace_back(std::stod(part));
	  ++group_id;
	}
      }
    }
  } else {
    throw std::runtime_error("Could not read from the likelihoods file.");
  }
  return seamat::DenseMatrix<double>(likelihoods);
}

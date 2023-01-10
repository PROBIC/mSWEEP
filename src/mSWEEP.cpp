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

void ReadLikelihoodFromFile(const std::string &likelihood_path, const Reference &reference, std::ostream &log, const std::unique_ptr<Sample> &sample) {
  log << "  reading likelihoods from file" << '\n';
  if (reference.get_n_groupings() > 1) {
    throw std::runtime_error("Using more than one grouping with --read-likelihood is not yet implemented.");
  }
  cxxio::In likelihoods(likelihood_path);
  sample->read_likelihood(reference.get_grouping(0), likelihoods.stream());
}

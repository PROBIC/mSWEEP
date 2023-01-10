#ifndef MSWEEP_MSWEEP_HPP
#define MSWEEP_MSWEEP_HPP

#include <vector>
#include <memory>
#include <cstddef>
#include <string>

#include "Sample.hpp"
#include "Reference.hpp"
#include "Grouping.hpp"

void ReadGroupIndicators(const std::string &indicators_path, Reference *reference);

void ReadPseudoalignments(const std::vector<std::string> &alignment_paths,
			  const std::string &themisto_merge_mode,
			  const Reference &reference,
			  std::unique_ptr<Sample> &sample);

void ReadLikelihoodFromFile(const std::string &likelihood_path, const Reference &reference, std::ostream &log, const std::unique_ptr<Sample> &sample);

#endif

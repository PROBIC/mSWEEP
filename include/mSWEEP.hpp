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
			  const bool compact_alignments,
			  const Reference &reference,
			  std::unique_ptr<Sample> &sample);

void ReadLikelihoodFromFile(const std::string &likelihood_path, const Reference &reference, std::ostream &log, const std::unique_ptr<Sample> &sample);

void ConstructLikelihood(const double tol,
			 const double frac_mu,
			 const Grouping &grouping,
			 const std::vector<uint32_t> &group_indicators,
			 const std::unique_ptr<Sample> &sample);

#endif

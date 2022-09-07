#ifndef MSWEEP_MSWEEP_HPP
#define MSWEEP_MSWEEP_HPP

#include <vector>
#include <memory>
#include <cstddef>

#include "parse_arguments.hpp"
#include "Sample.hpp"
#include "Reference.hpp"
#include "Grouping.hpp"

void ReadGroupIndicators(const Arguments &args, Reference *reference);

void ReadPseudoalignments(const Arguments &args,
			  const Reference &reference,
			  std::unique_ptr<Sample> &sample);

void ReadLikelihoodFromFile(const Arguments &args, const Reference &reference, std::ostream &log, const std::unique_ptr<Sample> &sample);

void ConstructLikelihood(const Arguments &args, const Grouping &grouping,
			 const std::vector<uint32_t> &group_indicators,
			 const std::unique_ptr<Sample> &sample, bool free_ec_counts);

void WriteResults(const Arguments &args, const std::unique_ptr<Sample> &sample,
		  const Grouping &grouping, const uint16_t n_groupings,
		  const uint16_t current_grouping);

#endif

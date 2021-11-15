#ifndef MSWEEP_MSWEEP_HPP
#define MSWEEP_MSWEEP_HPP

#include <vector>
#include <memory>
#include <cstddef>

#include "parse_arguments.hpp"
#include "Sample.hpp"
#include "Reference.hpp"
#include "Grouping.hpp"

void ReadInput(const Arguments &args,
	       std::vector<std::unique_ptr<Sample>> *samples,
	       std::ostream &log, Reference *reference);

void ConstructLikelihood(const Arguments &args, const Grouping &grouping,
			 const std::vector<uint32_t> &group_indicators,
			 const std::unique_ptr<Sample> &sample, bool free_ec_counts);

void WriteResults(const Arguments &args, const std::unique_ptr<Sample> &sample,
		  const Grouping &grouping, const uint16_t n_groupings,
		  const uint16_t current_grouping);

#endif

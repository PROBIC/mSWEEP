#ifndef MSWEEP_LIKELIHOOD_HPP
#define MSWEEP_LIKELIHOOD_HPP

#include "Grouping.hpp"
#include "Sample.hpp"

void likelihood_array_mat(const Grouping &grouping, const std::vector<uint32_t> &group_indicators, const double bb_constants[2], Sample &sample);

#endif

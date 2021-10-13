#ifndef MSWEEP_LIKELIHOOD_HPP
#define MSWEEP_LIKELIHOOD_HPP

#include "matrix.hpp"
#include "Reference.hpp"

void precalc_lls(const Grouping &grouping, const double bb_constants[2], Matrix<double> *ll_mat);

#endif

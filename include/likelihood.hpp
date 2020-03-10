#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include "matrix.hpp"
#include "Sample.hpp"
#include "Reference.hpp"

void precalc_lls(const Grouping &grouping, Matrix<double> *ll_mat);

#endif

#ifndef MSWEEP_LIKELIHOOD_HPP
#define MSWEEP_LIKELIHOOD_HPP

#include "Grouping.hpp"
#include "Sample.hpp"

seamat::DenseMatrix<double> likelihood_array_mat(const telescope::GroupedAlignment &pseudos, const Grouping &grouping, const double tol, const double frac_mu);

#endif

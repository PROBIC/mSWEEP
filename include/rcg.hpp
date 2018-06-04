#ifndef RCG_H
#define RCG_H

#include <vector>

#include "matrix.hpp"
#include "Sample.hpp"

Matrix<double> rcg_optl_mat(const Matrix<double> &logl, const Sample &sample, const std::vector<double> &alpha0, const double &tol, unsigned maxiters);

#endif

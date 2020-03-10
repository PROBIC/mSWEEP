#ifndef MSWEEP_RCG_HPP
#define MSWEEP_RCG_HPP

#include <vector>

#include "matrix.hpp"
#include "Sample.hpp"

Matrix<double> rcg_optl_mat(const Matrix<double> &logl, const Sample &sample, const std::vector<double> &alpha0, const double &tol, uint16_t maxiters);

#endif

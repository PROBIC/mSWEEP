#ifndef RCG_H
#define RCG_H

#include <vector>

#include "matrix.hpp"
#include "Sample.hpp"

Matrix<double> rcg_optl_mat(const Matrix<double> &logl, const Sample &sample, const std::vector<double> &alpha0, const double &tol, unsigned maxiters);
Matrix<double> rcg_optl_mat(const Matrix<double> &logl, const long unsigned &total_counts, const std::vector<double> &ec_counts, const std::vector<double> &alpha0, const double &tol, unsigned maxiters);

#endif

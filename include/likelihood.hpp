#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include "matrix.hpp"
#include "Sample.hpp"
#include "Reference.hpp"

Matrix<double> likelihood_array_mat(const Sample &sample, const Grouping &grouping);

#endif

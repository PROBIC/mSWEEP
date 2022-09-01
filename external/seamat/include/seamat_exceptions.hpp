// seamat: templatized matrix library
// https://github.com/tmaklin/seamat
//
// Copyright (C) 2021 Tommi Mäklin (tommi@maklin.fi)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef SEAMAT_SEAMAT_EXCEPTIONS_HPP
#define SEAMAT_SEAMAT_EXCEPTIONS_HPP

#define SEAMAT_CHECK_BOUNDS 1

#if defined(SEAMAT_CHECK_BOUNDS) && (SEAMAT_CHECK_BOUNDS) == 1
// Only verify dimensions if bounds checking is enabled

#include <exception>

#include "Matrix.hpp"

namespace seamat {
    template <typename T>
    const std::string DimensionsToString(const Matrix<T> &mat) {
	std::string msg(std::to_string(mat.get_rows()) + "x" + std::to_string(mat.get_cols()));
	return msg;
    }

    template <typename T>
    void MatrixSizesAreEqual(const Matrix<T> &lhs, const Matrix<T> &rhs) {
	if (lhs.get_rows() != rhs.get_rows() || lhs.get_cols() != rhs.get_cols()) {
	    std::string msg;
	    msg += "Matrix dimensions do not match:\n\t";
	    msg += "lhs: " + DimensionsToString(lhs) + "\n\t";
	    msg += "rhs: " + DimensionsToString(rhs) + "\n";
	    throw std::domain_error(msg);
	}
    }

    template <typename T>
    void MatricesCanBeMultiplied(const Matrix<T> &lhs, const Matrix<T> &rhs) {
	if (lhs.get_cols() != rhs.get_rows()) {
	    std::string msg;
	    msg += "Matrix dimensions do not allow multiplication:\n\t";
	    msg += "lhs: " + DimensionsToString(lhs) + "\n\t";
	    msg += "rhs: " + DimensionsToString(rhs) + "\n";
	    throw std::domain_error(msg);
	}
    }
    template <typename T, typename U>
    void MatrixCanBeMultipliedWithVector(const Matrix<T> &lhs, const std::vector<U> &rhs) {
	if (lhs.get_cols() != rhs.size()) {
	    std::string msg;
	    msg += "Matrix and vector dimensions do not allow multiplication:\n\t";
	    msg += DimensionsToString(lhs) + "\n\t";
	    msg += "rhs: " + std::to_string(rhs.size()) + "x1";
	    throw std::domain_error(msg);
	}
    }
}
#endif

#endif

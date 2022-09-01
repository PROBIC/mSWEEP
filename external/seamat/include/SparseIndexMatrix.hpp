// seamat: templatized matrix library
// https://github.com/tmaklin/seamat
//
// Copyright (C) 2021 Tommi Mäklin (tommi@maklin.fi)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef SEAMAT_SPARSE_INDEX_MATRIX_CPP
#define SEAMAT_SPARSE_INDEX_MATRIX_CPP

#include <cmath>
#include <stdexcept>

#include "Matrix.hpp"
#include "SparseIntegerTypeMatrix.hpp"
#include "openmp_config.hpp"
#include "math_util.hpp"

namespace seamat {
template <typename T, typename U> class SparseIndexMatrix : public Matrix<T> {
private:
    size_t n_rows_vals;
    size_t n_cols_vals;
    // Note: assumes that both vals and indices are sparse
    // For typical use the vals matrix is small enough that sparse/dense *probably* does not matter.
    // TODO: implement option to use dense vals.
    SparseMatrix<T> vals;
    SparseIntegerTypeMatrix<U> indices;

public:
    SparseIndexMatrix() = default;
    ~SparseIndexMatrix() = default;

    // Initialize from vals and indices
    SparseIndexMatrix(const Matrix<T> &_vals, const Matrix<U> &_indices, const bool store_as_sparse);

    // Access individual elements
    T& operator()(size_t row, size_t col) override;
    const T& operator()(size_t row, size_t col) const override;

    // Mathematical operators
    // Matrix-matrix in-place summation and subtraction
    SparseIndexMatrix<T, U>& operator+=(const Matrix<T>& rhs) override;
    SparseIndexMatrix<T, U>& operator-=(const Matrix<T>& rhs) override;

    // In-place right multiplication
    SparseIndexMatrix<T, U>& operator*=(const Matrix<T>& rhs) override;
    // In-place left multiplication
    SparseIndexMatrix<T, U>& operator%=(const Matrix<T>& rhs) override;

    // Matrix-scalar, in-place
    SparseIndexMatrix<T, U>& operator+=(const T& rhs) override;
    SparseIndexMatrix<T, U>& operator-=(const T& rhs) override;
    SparseIndexMatrix<T, U>& operator*=(const T& rhs) override;
    SparseIndexMatrix<T, U>& operator/=(const T& rhs) override;

};

// Access individual elements
template <typename T, typename U>
T& SparseIndexMatrix<T,U>::operator()(size_t row, size_t col) {
    size_t out_col = (*this->indices)(row, col);
    return (*this->vals)(row, out_col);
}

// Access individual elements (const)
template <typename T, typename U>
const T& SparseIndexMatrix<T,U>::operator()(size_t row, size_t col) const {
    size_t out_col = (*this->indices)(row, col);
    return (*this->vals)(row, out_col);
}

// Initialize from vals and indices
template<typename T, typename U>
SparseIndexMatrix<T, U>::SparseIndexMatrix(const Matrix<T> &_vals, const Matrix<U> &_indices, const bool store_as_sparse) {
    this->resize_rows(_indices.get_rows());
    this->resize_cols(_indices.get_cols());

    this->n_rows_vals = _vals.get_rows();
    this->n_cols_vals = _vals.get_cols();
    if (store_as_sparse) {
	this->indices.reset(new SparseMatrix<U>(_indices, 0));
	this->vals.reset(new SparseMatrix<T>(_vals, 0.0));
    } else {
	this->indices.reset(new DenseMatrix<U>(_indices));
	this->vals.reset(new DenseMatrix<T>(_vals));
    }
}

// TODO implement index matrix operators

// In-place matrix-matrix addition
template<typename T, typename U>
SparseIndexMatrix<T, U>& SparseIndexMatrix<T, U>::operator+=(const Matrix<T>& rhs) {
    // seamat::SparseIndexMatrix<T, U>::operator+=
    //
    // Add values from rhs to the calling matrix in-place.
    //
    //   Input:
    //     `rhs`: Matrix to add, must have the same dimensions as the caller.
    //
    throw std::runtime_error("Index matrix operators have not been implemented.");
    return *this;
}

// In-place matrix-matrix subtraction
template<typename T, typename U>
SparseIndexMatrix<T, U>& SparseIndexMatrix<T, U>::operator-=(const Matrix<T>& rhs) {
    // seamat::SparseIndexMatrix<T, U>::operator-=
    //
    // Subtract values of rhs from the calling matrix in-place.
    //
    //   Input:
    //     `rhs`: Matrix to subtract, must have the same dimensions as the caller.
    //
    throw std::runtime_error("Index matrix operators have not been implemented.");

    return *this;
}

// In-place right multiplication
template<typename T, typename U>
SparseIndexMatrix<T, U>& SparseIndexMatrix<T, U>::operator*=(const Matrix<T>& rhs) {
    // seamat::SparseIndexMatrix<T, U>::operator*=
    //
    // Matrix right-multiplication of the caller with rhs in-place.
    //
    //   Input:
    //     `rhs`: Matrix to right multiply with,
    //            must have the same number of rows as the caller has columns.
    //
    throw std::runtime_error("Index matrix operators have not been implemented.");
    return *this;
}

// In-place left multiplication
template<typename T, typename U>
SparseIndexMatrix<T, U>& SparseIndexMatrix<T, U>::operator%=(const Matrix<T>& lhs) {
    // seamat::SparseIndexMatrix<T, U>::operator%=
    //
    // Matrix left-multiplication of the caller with lhs in-place.
    //
    //   Input:
    //     `lhs`: Matrix to left multiply with,
    //            must have the same number of columns as the caller has rows.
    //
    throw std::runtime_error("Index matrix operators have not been implemented.");
    return *this;
}

// In-place matrix-scalar addition
template<typename T, typename U>
SparseIndexMatrix<T, U>& SparseIndexMatrix<T, U>::operator+=(const T& scalar) {
    // seamat::SparseIndexMatrix<T, U>::operator+=
    //
    // In-place addition of a scalar to caller.
    //
    //   Input:
    //     `scalar`: Scalar value to add to all caller values.
    //
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->n_rows_vals; ++i) {
	for (size_t j = 0; j < this->n_cols_vals; ++j) {
	    this->vals->operator()(i, j) += scalar;
	}
    }

    return *this;
}

// In-place matrix-scalar subtraction
template<typename T, typename U>
SparseIndexMatrix<T, U>& SparseIndexMatrix<T, U>::operator-=(const T& scalar) {
    // seamat::SparseIndexMatrix<T, U>::operator-=
    //
    // In-place subtraction of a scalar from the caller.
    //
    //   Input:
    //     `scalar`: Scalar value to subtract from all caller values.
    //
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->n_rows_vals; ++i) {
	for (size_t j = 0; j < this->n_cols_vals; ++j) {
	    this->vals->operator()(i, j) -= scalar;
	}
    }

    return *this;
}

// In-place matrix-scalar multiplication
template<typename T, typename U>
SparseIndexMatrix<T, U>& SparseIndexMatrix<T, U>::operator*=(const T& scalar) {
    // seamat::SparseIndexMatrix<T, U>::operator*=
    //
    // In-place multiplication of the caller with a scalar.
    //
    //   Input:
    //     `scalar`: Scalar value to multiply all caller values with.
    //
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->n_rows_vals; ++i) {
	for (size_t j = 0; j < this->n_cols_vals; ++j) {
	    this->vals->operator()(i, j) *= scalar;
	}
    }

    return *this;
}

// In-place matrix-scalar division
template<typename T, typename U>
SparseIndexMatrix<T, U>& SparseIndexMatrix<T, U>::operator/=(const T& scalar) {
    // seamat::SparseIndexMatrix<T, U>::operator/=
    //
    // In-place division of the caller with a scalar.
    //
    //   Input:
    //     `scalar`: Scalar value to divide all caller values with.
    //
    if (nearly_equal<T>(scalar, (T)0))
	throw std::runtime_error("Math error: attempt to divide by zero.");

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->n_rows_vals; ++i) {
	for (size_t j = 0; j < this->n_cols_vals; ++j) {
	    this->vals->operator()(i, j) /= scalar;
	}
    }

    return *this;
}
}

#endif

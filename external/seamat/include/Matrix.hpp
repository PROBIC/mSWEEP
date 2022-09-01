// seamat: templatized matrix library
// https://github.com/tmaklin/seamat
//
// Copyright (C) 2021 Tommi Mäklin (tommi@maklin.fi)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef SEAMAT_MATRIX_HPP
#define SEAMAT_MATRIX_HPP

#include <vector>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>

#include "bm.h"
#include "bmsparsevec.h"

// Basic matrix structure and operations
// Implementation was done following the instructions at
// https://www.quantstart.com/articles/Matrix-Classes-in-C-The-Header-File
//
// **None of the operations validate the matrix sizes**
//
// This file provides implementations for the functions that are not
// parallellized.  Implementations for the functions that are
// parallellized should be provided in the .cpp file that is included
// at the very end of this file.

namespace seamat {
// Forward declare implementations of abstract class Matrix
template <typename T> class DenseMatrix;
template <typename T> class SparseMatrix;

template <typename T> class Matrix {
private:
    size_t rows;
    size_t cols;

protected:
    // Derived classes can use these to resize the base class
    void resize_rows(const size_t new_rows) { this->rows = new_rows; };
    void resize_cols(const size_t new_cols) { this->cols = new_cols; };

public:
    //////
    // Pure virtuals - must override in derived classes.
    ////
    // Access individual elements
    virtual T& operator()(size_t row, size_t col) =0;
    virtual const T& operator()(size_t row, size_t col) const =0;

    // Mathematical operators
    // Matrix-matrix in-place summation and subtraction
    virtual Matrix<T>& operator+=(const Matrix<T>& rhs) =0;
    virtual Matrix<T>& operator-=(const Matrix<T>& rhs) =0;

    // In-place right multiplication
    virtual Matrix<T>& operator*=(const Matrix<T>& rhs) =0;
    // In-place left multiplication
    virtual Matrix<T>& operator%=(const Matrix<T>& rhs) =0;

    // Matrix-scalar, in-place
    virtual Matrix<T>& operator+=(const T& rhs) =0;
    virtual Matrix<T>& operator-=(const T& rhs) =0;
    virtual Matrix<T>& operator*=(const T& rhs) =0;
    virtual Matrix<T>& operator/=(const T& rhs) =0;

    //////
    // Implemented in Matrix.cpp - these functions only use the pure
    // virtual functions and the private member variables, and are
    // generic in the sense that they work with any input class
    // derived from Matrix<T>.
    ////
    // Mathematical operators
    // Matrix-matrix summation and subtraction
    DenseMatrix<T>& operator+(const Matrix<T>& rhs) const;
    DenseMatrix<T>& operator-(const Matrix<T>& rhs) const;

    // Matrix-matrix right multiplication
    DenseMatrix<T> operator*(const Matrix<T>& rhs) const;
    // Matrix-matrix left multiplication
    DenseMatrix<T> operator%(const Matrix<T>& rhs) const;

    // Matrix-scalar
    DenseMatrix<T> operator+(const T& rhs);
    DenseMatrix<T> operator-(const T& rhs);
    DenseMatrix<T> operator*(const T& rhs);
    DenseMatrix<T> operator/(const T& rhs);

    // Matrix-matrix comparisons
    bool operator==(const Matrix<T>& rhs) const;

    // Matrix-vector multiplication
    template <typename U>
    std::vector<T> operator*(const std::vector<U>& rhs) const;
    template <typename U>
    std::vector<U> operator*(const std::vector<U>& rhs) const;

    // Matrix-vector right multiplication, store result in arg
    template <typename U, typename V>
    void right_multiply(const std::vector<U>& rhs, std::vector<V> &result) const;
    template <typename U, typename V>
    void logspace_right_multiply(const std::vector<U>& rhs, std::vector<V>& result) const;

    // LogSumExp a column
    template <typename V>
    T log_sum_exp_col(const size_t col_id) const;

    // Generic transpose
    DenseMatrix<T> transpose() const;

    // Get the number of rows in the matrix
    size_t get_rows() const { return this->rows; }
    // Get the number of columns of the matrix
    size_t get_cols() const { return this->cols; }

    // Turn an arbitrary matrix into a dense matrix
    DenseMatrix<T> densify() const { return DenseMatrix<T>(this); }
    // Turn an arbitrary matrix into a column sparse matrix
    SparseMatrix<T> sparsify(const T &zero_val) const { return SparseMatrix<T>(this, zero_val); }

};
// Matrix-matrix addition
template<typename T>
DenseMatrix<T>& Matrix<T>::operator+(const Matrix<T>& rhs) const {
    // seamat::Matrix<T>::operator+
    //
    // Creates a new DenseMatrix<T> that contains the result of
    // summing two instances of Matrix<T>.
    //
    //   Input:
    //     `rhs`: Matrix to add,
    //            must have the same dimensions as the caller.
    //
    //   Output:
    //   `result`: A new matrix containing the sum of the caller and `rhs`.
    //
#if defined(SEAMAT_CHECK_BOUNDS) && (SEAMAT_CHECK_BOUNDS) == 1
    try {
	MatrixSizesAreEqual(*this, rhs);
    } catch (const std::exception &e) {
	throw e;
    }
#endif

    DenseMatrix<T> result(this->get_rows(), this->get_cols(), (T)0);

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    result(i, j) = this->operator()(i, j) + rhs(i, j);
	}
    }

    return result;
}

// Matrix-matrix subtraction
template<typename T>
DenseMatrix<T>& Matrix<T>::operator-(const Matrix<T>& rhs) const {
    // seamat::Matrix<T>::operator-
    //
    // Creates a new DenseMatrix<T> that contains the result of
    // subtracting rhs from the caller.
    //
    //   Input:
    //     `rhs`: Matrix to subtract,
    //            must have the same dimensions as the caller.
    //
    //   Output:
    //   `result`: A new matrix containing the result.
    //
#if defined(SEAMAT_CHECK_BOUNDS) && (SEAMAT_CHECK_BOUNDS) == 1
    try {
	MatrixSizesAreEqual(*this, rhs);
    } catch (const std::exception &e) {
	throw e;
    }
#endif

    DenseMatrix<T> result(this->get_rows(), this->get_cols(), (T)0);

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    result(i, j) = this->operator()(i, j) - rhs(i, j);
	}
    }

    return result;
}

// Matrix-matrix right multiplication
template<typename T>
DenseMatrix<T> Matrix<T>::operator*(const Matrix<T>& rhs) const {
    // seamat::Matrix<T>::operator*
    //
    // Creates a new DenseMatrix<T> that contains the result of
    // right multiplying the caller with rhs.
    //
    //   Input:
    //     `rhs`: Matrix to right multiply with,
    //            must have the same number of rows as the caller has columns.
    //
    //   Output:
    //   `result`: A new matrix containing the result.
    //
#if defined(SEAMAT_CHECK_BOUNDS) && (SEAMAT_CHECK_BOUNDS) == 1
    try {
	MatricesCanBeMultiplied(*this, rhs);
    } catch (const std::exception &e) {
	throw e;
    }
#endif

    DenseMatrix<T> result(this->get_rows(), rhs.get_cols(), (T)0);

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    for (size_t k = 0; k < rhs.get_cols(); k++) {
		result(i, k) += this->operator()(i, j) * rhs(j, k);
	    }
	}
    }

    return result;
}

// Matrix-matrix left multiplication
template<typename T>
DenseMatrix<T> Matrix<T>::operator%(const Matrix<T>& lhs) const {
    // seamat::Matrix<T>::operator%
    //
    // Creates a new DenseMatrix<T> that contains the result of
    // left multiplying the caller with lhs.
    //
    //   Input:
    //     `lhs`: Matrix to right multiply with,
    //            must have the same number of columns as the caller has rows.
    //
    //   Output:
    //   `result`: A new matrix containing the result.
    //
    DenseMatrix<T> &result = lhs * (*this); // operator* checks bounds
    return result;
}

// Matrix-scalar addition
template<typename T>
DenseMatrix<T> Matrix<T>::operator+(const T& scalar) {
    // seamat::Matrix<T>::operator+
    //
    // Addition of a scalar to caller.
    //
    //   Input:
    //     `scalar`: Scalar value to add to all caller values.
    //
    //   Output:
    //     `result`: Result of the addition.
    //
    DenseMatrix<T> result(this->get_rows(), this->get_cols(), scalar);
    result += *this;

    return result;
}

// Matrix-scalar subtraction
template<typename T>
DenseMatrix<T> Matrix<T>::operator-(const T& scalar) {
    // seamat::Matrix<T>::operator-
    //
    // Subtraction of a scalar from the caller.
    //
    //   Input:
    //     `scalar`: Scalar value to subtract from all caller values.
    //
    //   Output:
    //     `result`: Result of the subtraction.
    //
    DenseMatrix<T> result(this->get_rows(), this->get_cols(), -scalar);
    result += *this;

    return result;
}

// Matrix-scalar multiplication
template<typename T>
DenseMatrix<T> Matrix<T>::operator*(const T& scalar) {
    // seamat::Matrix<T>::operator*
    //
    // Multiplication of the caller with a scalar.
    //
    //   Input:
    //     `scalar`: Scalar value to multiply all caller values with.
    //
    //   Output:
    //     `result`: Result of the multiplication.
    //
    DenseMatrix<T> result(this->get_rows(), this->get_cols(), (T)0);
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    result(i, j) = this->operator()(i, j) * scalar;
	}
    }

    return result;
}

// Matrix-scalar division
template<typename T>
DenseMatrix<T> Matrix<T>::operator/(const T& scalar) {
    // seamat::Matrix<T>::operator/
    //
    // Division of the caller with a scalar.
    //
    //   Input:
    //     `scalar`: Scalar value to divide all caller values with.
    //
    //   Output:
    //     `result`: Result of the division.
    //
    DenseMatrix<T> result(this->get_rows(), this->get_cols(), (T)0);
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    result(i, j) = this->operator()(i, j) / scalar;
	}
    }

    return result;
}

// Matrix-matrix comparison
template<typename T>
bool Matrix<T>::operator==(const Matrix<T>& rhs) const {
    // seamat::Matrix<T>::operator==
    //
    // Check if two matrices contain exactly the same values.
    // Note: two matrices are not necessarily identical objects
    // even if they contain the same values.
    //
    //   Input:
    //     `rhs`: Matrix to check values with.
    //
    // TODO: implement a specialization for floating point numbers using std::nextafter.
    //
#if defined(SEAMAT_CHECK_BOUNDS) && (SEAMAT_CHECK_BOUNDS) == 1
    try {
	MatrixSizesAreEqual(*this, rhs);
    } catch (const std::exception &e) {
	throw e;
    }
#endif
    bool all_equal = true;
#pragma omp parallel for schedule(static) reduction (&:all_equal)
    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    all_equal &= (this->operator()(i, j) == rhs(i, j));
	}
    }

    return all_equal;
}

// Matrix-vector right multiplication, return same type as caller
template<typename T>
template <typename U>
std::vector<T> Matrix<T>::operator*(const std::vector<U>& rhs) const {
    // seamat::Matrix<T>::operator*
    //
    // Right-multiply the caller with vector rhs.
    //
    //   Input:
    //     `rhs`: Vector to right multiply with. Must have the same number
    //            of elements as the caller has columns.
    //
    //   Output:
    //     `result`: Vector containing the result of the right multiplication.
    //
#if defined(SEAMAT_CHECK_BOUNDS) && (SEAMAT_CHECK_BOUNDS) == 1
    try {
	MatrixCanBeMultipliedWithVector(*this, rhs);
    } catch (const std::exception &e) {
	throw e;
    }
#endif

    std::vector<T> result(this->get_rows(), 0.0);

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    result[i] += this->operator()(i, j) * rhs[j];
	}
    }

    return result;
}

// Matrix-vector right multiplication, return same type as rhs
template<typename T>
template <typename U>
std::vector<U> Matrix<T>::operator*(const std::vector<U>& rhs) const {
    // seamat::Matrix<T>::operator*
    //
    // Right-multiply the caller with vector rhs.
    //
    //   Input:
    //     `rhs`: Vector to right multiply with. Must have the same number
    //            of elements as the caller has columns.
    //
    //   Output:
    //     `result`: Vector containing the result of the right multiplication.
    //
#if defined(SEAMAT_CHECK_BOUNDS) && (SEAMAT_CHECK_BOUNDS) == 1
    try {
	MatrixCanBeMultipliedWithVector(*this, rhs);
    } catch (const std::exception &e) {
	throw e;
    }
#endif

    std::vector<U> result(this->get_rows(), 0.0);

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    result[i] += this->operator()(i, j) * rhs[j];
	}
    }

    return result;
}

// Matrix-vector right multiplication, store result in argument
template<typename T>
template <typename U, typename V>
void Matrix<T>::right_multiply(const std::vector<U>& rhs, std::vector<V> &result) const {
    // seamat::Matrix<T>::right_multiply
    //
    // Right-multiply the caller with vector rhs and store result in argument.
    //
    //   Input:
    //     `rhs`: Vector to right multiply with. Must have the same number
    //            of elements as the caller has columns.
    //     `result:` Vector to store the result in. Must have the same number
    //               of elements as the caller has columns.
    //
#if defined(SEAMAT_CHECK_BOUNDS) && (SEAMAT_CHECK_BOUNDS) == 1
    try {
	MatrixCanBeMultipliedWithVector(*this, rhs);
	MatrixCanBeMultipliedWithVector(*this, result);
    } catch (const std::exception &e) {
	throw e;
    }
#endif

    std::fill(result.begin(), result.end(), (V)0);
#pragma omp parallel for schedule(static) reduction(vec_double_plus:result)
    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    result[i] += this->operator()(i, j) * rhs[j];
	}
    }
}

// log-space Matrix-vector right multiplication, store result in arg
template<typename T>
template <typename U, typename V>
void Matrix<T>::logspace_right_multiply(const std::vector<U>& rhs, std::vector<V>& result) const {
    // seamat::Matrix<T>::logspace_right_multiply
    //
    // Right-multiply the caller with vector rhs in logspace and store result in argument.
    //
    //   Input:
    //     `rhs`: Vector to right multiply with. Must have the same number
    //            of elements as the caller has columns.
    //     `result:` Vector to store the result in. Must have the same number
    //               of elements as the caller has columns.
    //
#if defined(SEAMAT_CHECK_BOUNDS) && (SEAMAT_CHECK_BOUNDS) == 1
    try {
	MatrixCanBeMultipliedWithVector(*this, rhs);
	MatrixCanBeMultipliedWithVector(*this, result);
    } catch (const std::exception &e) {
	throw e;
    }
#endif

    std::fill(result.begin(), result.end(), (V)0);
#pragma omp parallel for schedule(static) reduction(vec_double_plus:result)
    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    result[i] += std::exp(this->operator()(i, j) + rhs[j]);
	}
    }
}

// LogSumExp a column
template<typename T>
template <typename V> // desired accuracy (float, double, long double, etc.)
T Matrix<T>::log_sum_exp_col(const size_t col_id) const {
    // seamat::Matrix<T>::log_sum_exp_col
    //
    // Calculates the log of the sum of exponentials of the values in column col_id.
    //
    //   Input:
    //     `col_id`: which column to LogSumExp.
    //
    //   Output:
    //     `max_elem + std::log(sum)`: Log of the sum of exponentials of values in column col_id.
    //
    // Note: this function accesses the elements inefficiently since
    //       the matrices are implement row-wise. This means that it
    //       shouldn't be parallellised here, however, the caller can
    //       parallellize logsumexping multiple cols.
    //
#if defined(SEAMAT_CHECK_BOUNDS) && (SEAMAT_CHECK_BOUNDS) == 1
    if (col_id >= this->get_cols())
	throw std::domain_error("Column " + std::to_string(col_id) + " is out of bounds:\n\t Matrix dimensions are: " + DimensionsToString(*this) + ".\n");
#endif

    T max_elem = 0;
    V exp_sum = 0;
    for (size_t i = 0; i < this->get_rows(); ++i) {
	max_elem = (this->operator()(i, col_id) > max_elem ? this->operator()(i, col_id) : max_elem);
    }

    for (size_t i = 0; i < this->get_rows(); ++i) {
	exp_sum += std::exp((V)this->operator()(i, col_id) - (V)max_elem);
    }
    return (V)max_elem + std::log(exp_sum);
}

// Generic transpose
template<typename T>
DenseMatrix<T> Matrix<T>::transpose() const {
    // seamat::Matrix<T>::transpose
    //
    // Transpose a matrix.
    //
    //   Output:
    //     `result`: Transpose of the caller as a dense matrix.
    //
    DenseMatrix<T> result(this->get_rows(), this->get_cols(), (T)0);

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    result(i, j) = this->operator()(j, i);
	}
    }

    return result;
}
}

#include "DenseMatrix.hpp"
#include "IndexMatrix.hpp"
#include "SparseIndexMatrix.hpp"
#include "SparseIntegerTypeMatrix.hpp"
#include "SparseMatrix.hpp"

#endif

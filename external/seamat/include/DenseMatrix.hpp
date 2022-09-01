// seamat: templatized matrix library
// https://github.com/tmaklin/seamat
//
// Copyright (C) 2021 Tommi Mäklin (tommi@maklin.fi)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef SEAMAT_DENSE_MATRIX_CPP
#define SEAMAT_DENSE_MATRIX_CPP

#include <cmath>

#include "Matrix.hpp"
#include "openmp_config.hpp"

namespace seamat {
template <typename T> class DenseMatrix : public Matrix<T> {
private:
    std::vector<T> mat;

public:
    DenseMatrix() = default;
    ~DenseMatrix() = default;
    // Parameter constructor
    DenseMatrix(size_t _rows, size_t _cols, const T& _initial);
    // Copy constructor from contiguous 2D vector
    DenseMatrix(const std::vector<T> &rhs, const size_t _rows, const size_t _cols);
    // Copy constructor from 2D vector
    DenseMatrix(const std::vector<std::vector<T>> &rhs);
    // Copy constructor from another matrix
    DenseMatrix(const Matrix<T> &rhs);

    // Assignment operator
    DenseMatrix<T>& operator=(const Matrix<T>& rhs);

   // Resize a matrix
    void resize(const size_t new_rows, const size_t new_cols, const T initial);

    // Access individual elements
    T& operator()(size_t row, size_t col) override;
    const T& operator()(size_t row, size_t col) const override;

    // Mathematical operators
    // Matrix-matrix in-place summation and subtraction
    DenseMatrix<T>& operator+=(const Matrix<T>& rhs) override;
    DenseMatrix<T>& operator-=(const Matrix<T>& rhs) override;

    // In-place right multiplication
    DenseMatrix<T>& operator*=(const Matrix<T>& rhs) override;
    // In-place left multiplication
    DenseMatrix<T>& operator%=(const Matrix<T>& rhs) override;

    // Matrix-scalar, in-place
    DenseMatrix<T>& operator+=(const T& rhs);
    DenseMatrix<T>& operator-=(const T& rhs);
    DenseMatrix<T>& operator*=(const T& rhs);
    DenseMatrix<T>& operator/=(const T& rhs);

    // Fill a matrix with the sum of two matrices in-place
    template <typename V, typename U>
    void sum_fill(const Matrix<V>& rhs1, const Matrix<U>& rhs2);

};

// Parameter Constructor
template<typename T>
DenseMatrix<T>::DenseMatrix(size_t _rows, size_t _cols, const T& _initial) {
    mat.resize(_rows*_cols, _initial);
    this->resize_rows(_rows);
    this->resize_cols(_cols);
}

// Copy constructor from contiguous 2D vector
template <typename T>
DenseMatrix<T>::DenseMatrix(const std::vector<T> &rhs, const size_t _rows, const size_t _cols) {
    mat = rhs;
    this->resize_rows(_rows);
    this->resize_cols(_cols);
}

// Copy constructor from 2D vector
template<typename T>
DenseMatrix<T>::DenseMatrix(const std::vector<std::vector<T>> &rhs) {
    this->resize_rows(rhs.size());
    this->resize_cols(rhs.at(0).size());
    mat.resize(this->get_rows()*this->get_cols());
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    this->operator()(i, j) = rhs[i][j];
	}
    }
}

// Copy constructor from another matrix
template<typename T>
DenseMatrix<T>::DenseMatrix(const Matrix<T> &rhs) {
    this->resize_rows(rhs.get_rows());
    this->resize_cols(rhs.get_cols());
    mat.resize(this->get_rows()*this->get_cols());
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    this->operator()(i, j) = rhs(i, j);
	}
    }
}

// Assignment Operator
template<typename T>
DenseMatrix<T>& DenseMatrix<T>::operator=(const Matrix<T>& rhs) {
    if (&rhs == this)
	return *this;

    size_t new_rows = rhs.get_rows();
    size_t new_cols = rhs.get_cols();
    if (new_rows != this->get_rows() || new_cols != this->get_cols()) {
	resize(new_rows, new_cols, (T)0);
    }
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < new_rows; i++) {
	for (size_t j = 0; j < new_cols; j++) {
	    this->operator()(i, j) = rhs(i, j);
	}
    }
    return *this;
}

// Resize a matrix
template<typename T>
void DenseMatrix<T>::resize(const size_t new_rows, const size_t new_cols, const T initial) {
    if (new_rows != this->get_rows() || new_cols != this->get_cols()) {
	mat.resize(new_rows*new_cols, initial);
	this->resize_rows(new_rows);
	this->resize_cols(new_cols);
    }
}

// Access individual elements
template <typename T>
T& DenseMatrix<T>::operator()(size_t row, size_t col) {
    return this->mat[row*this->get_cols() + col];
}

// Access individual elements (const)
template <typename T>
const T& DenseMatrix<T>::operator()(size_t row, size_t col) const {
    return this->mat[row*this->get_cols() + col];
}


// In-place matrix-matrix subtraction
template<typename T>
DenseMatrix<T>& DenseMatrix<T>::operator-=(const Matrix<T>& rhs) {
    // seamat::DenseMatrix<T>::operator-=
    //
    // Subtract values of rhs from the calling matrix in-place.
    //
    //   Input:
    //     `rhs`: Matrix to subtract, must have the same dimensions as the caller.
    //
#if defined(SEAMAT_CHECK_BOUNDS) && (SEAMAT_CHECK_BOUNDS) == 1
    try {
	MatrixSizesAreEqual(*this, rhs);
    } catch (const std::exception &e) {
	throw e;
    }
#endif

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->get_rows(); i++) {
	for (size_t j = 0; j < this->get_cols(); j++) {
	    this->operator()(i, j) -= rhs(i, j);
	}
    }

    return *this;
}

// In-place matrix-matrix addition
template<typename T>
DenseMatrix<T>& DenseMatrix<T>::operator+=(const Matrix<T>& rhs) {
    // seamat::DenseMatrix<T>::operator+=
    //
    // Add values from rhs to the calling matrix in-place.
    //
    //   Input:
    //     `rhs`: Matrix to add, must have the same dimensions as the caller.
    //
#if defined(SEAMAT_CHECK_BOUNDS) && (SEAMAT_CHECK_BOUNDS) == 1
    try {
	MatrixSizesAreEqual(*this, rhs);
    } catch (const std::exception &e) {
	throw e;
    }
#endif

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->get_rows(); i++) {
	for (size_t j = 0; j < this->get_cols(); j++) {
	    this->operator()(i, j) += rhs(i, j);
	}
    }

    return *this;
}

// In-place right multiplication
template<typename T>
DenseMatrix<T>& DenseMatrix<T>::operator*=(const Matrix<T>& rhs) {
    // seamat::DenseMatrix<T>::operator*=
    //
    // Matrix right-multiplication of the caller with rhs in-place.
    //
    //   Input:
    //     `rhs`: Matrix to right multiply with,
    //            must have the same number of rows as the caller has columns.
    //
    const DenseMatrix<T> &result = (*this) * rhs;
    (*this) = result;
    return *this;
}

// In-place left multiplication
template<typename T>
DenseMatrix<T>& DenseMatrix<T>::operator%=(const Matrix<T>& lhs) {
    // seamat::DenseMatrix<T>::operator%=
    //
    // Matrix left-multiplication of the caller with lhs in-place.
    //
    //   Input:
    //     `lhs`: Matrix to left multiply with,
    //            must have the same number of columns as the caller has rows.
    //
    const DenseMatrix<T> &result = lhs * (*this);
    (*this) = result;
    return *this;
}

// In-place matrix-scalar addition
template<typename T>
DenseMatrix<T>& DenseMatrix<T>::operator+=(const T& scalar) {
    // seamat::DenseMatrix<T>::operator+=
    //
    // In-place addition of a scalar to caller.
    //
    //   Input:
    //     `scalar`: Scalar value to add to all caller values.
    //
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->get_rows(); i++) {
	for (size_t j = 0; j < this->get_cols(); j++) {
	    this->operator()(i, j) += scalar;
	}
    }

    return *this;
}

// In-place matrix-scalar subtraction
template<typename T>
DenseMatrix<T>& DenseMatrix<T>::operator-=(const T& scalar) {
    // seamat::DenseMatrix<T>::operator-=
    //
    // In-place subtraction of a scalar from the caller.
    //
    //   Input:
    //     `scalar`: Scalar value to subtract from all caller values.
    //
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->get_rows(); i++) {
	for (size_t j=0; j < this->get_cols(); j++) {
	    this->operator()(i, j) -= scalar;
	}
    }

    return *this;
}

// In-place matrix-scalar multiplication
template<typename T>
DenseMatrix<T>& DenseMatrix<T>::operator*=(const T& scalar) {
    // seamat::DenseMatrix<T>::operator*=
    //
    // In-place multiplication of the caller with a scalar.
    //
    //   Input:
    //     `scalar`: Scalar value to multiply all caller values with.
    //
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    this->operator()(i, j) *= scalar;
	}
    }

    return *this;
}

// In-place matrix-scalar division
template<typename T>
DenseMatrix<T>& DenseMatrix<T>::operator/=(const T& scalar) {
    // seamat::DenseMatrix<T>::operator/=
    //
    // In-place division of the caller with a scalar.
    //
    //   Input:
    //     `scalar`: Scalar value to divide all caller values with.
    //
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    this->operator()(i, j) /= scalar;
	}
    }

    return *this;
}

// Fill matrix with sum of two matrices
template <typename T>
template <typename V, typename U>
void DenseMatrix<T>::sum_fill(const Matrix<V>& rhs1, const Matrix<U>& rhs2) {
#if defined(SEAMAT_CHECK_BOUNDS) && (SEAMAT_CHECK_BOUNDS) == 1
    try {
	MatrixSizesAreEqual(*this, rhs1);
	MatrixSizesAreEqual(*this, rhs2);
    } catch (const std::exception &e) {
	throw e;
    }
#endif

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    this->operator()(i, j) = (T)rhs1(i, j) + (T)rhs2(i, j);
	}
    }
}
}

#endif

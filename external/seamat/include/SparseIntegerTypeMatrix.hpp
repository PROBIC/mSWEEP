// seamat: templatized matrix library
// https://github.com/tmaklin/seamat
//
// Copyright (C) 2021 Tommi Mäklin (tommi@maklin.fi)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef SEAMAT_SPARSE_INTEGER_TYPE_MATRIX_CPP
#define SEAMAT_SPARSE_INTEGER_TYPE_MATRIX_CPP

#include <cmath>
#include <stdexcept>

#include "Matrix.hpp"
#include "openmp_config.hpp"
#include "math_util.hpp"

namespace seamat {
template <typename T> class SparseIntegerTypeMatrix : public Matrix<T> {
private:
    // Specialization of SparseIntegerTypeMatrix for integer types using BitMagic's sparse integer vectors
    //
    bm::sparse_vector<T, bm::bvector<>> vals;
    T zero_val;

public:
    SparseIntegerTypeMatrix() = default;
    ~SparseIntegerTypeMatrix() = default;
    // Parameter constructor
    SparseIntegerTypeMatrix(size_t _rows, size_t _cols, const T& _initial);
    // Initialize from a DenseMatrix
    SparseIntegerTypeMatrix(const Matrix<T> &_vals, const T& _zero_val);
    // Initialize from a 2D vector
    SparseIntegerTypeMatrix(const std::vector<std::vector<T>> &rhs, const T& _zero_val);
    // Copy constructor from contiguous 2D vector
    SparseIntegerTypeMatrix(const std::vector<T> &rhs, const size_t _rows, const size_t _cols, const T& _zero_val);

    // Access individual elements
    T& operator()(size_t row, size_t col) override;
    const T& operator()(size_t row, size_t col) const override;

    // Mathematical operators
    // Matrix-matrix in-place summation and subtraction
    SparseIntegerTypeMatrix<T>& operator+=(const Matrix<T>& rhs) override;
    SparseIntegerTypeMatrix<T>& operator-=(const Matrix<T>& rhs) override;

    // In-place right multiplication
    SparseIntegerTypeMatrix<T>& operator*=(const Matrix<T>& rhs) override;
    // In-place left multiplication
    SparseIntegerTypeMatrix<T>& operator%=(const Matrix<T>& rhs) override;

    // Matrix-scalar, in-place
    SparseIntegerTypeMatrix<T>& operator+=(const T& rhs) override;
    SparseIntegerTypeMatrix<T>& operator-=(const T& rhs) override;
    SparseIntegerTypeMatrix<T>& operator*=(const T& rhs) override;
    SparseIntegerTypeMatrix<T>& operator/=(const T& rhs) override;
};

// Parameter constructor
template<typename T>
SparseIntegerTypeMatrix<T>::SparseIntegerTypeMatrix(size_t _rows, size_t _cols, const T& _initial) {
    // Initializes a dense SparseIntegerTypeMatrix; use
    // SparseIntegerTypeMatrix<T>::remove_nonzeros to sparsify the matrix after
    // filling it.
    this->resize_rows(_rows);
    this->resize_cols(_cols);

    size_t n_elements = _rows*_cols;
    this->vals.resize(n_elements);
}
    
// Initialize from a DenseMatrix
template<typename T>
SparseIntegerTypeMatrix<T>::SparseIntegerTypeMatrix(const Matrix<T> &_vals, const T& _zero_val) {
    // TODO: get the values in _vals as an array and initialize the SparseMatrix from that?
    // According to the BitMagic sparsevector documentation this is the fastest way to initialize
    // a sparsevector type (calling bm::sparse_vector::set is supposedly slow).
    this->resize_rows(_vals.get_rows());
    this->resize_cols(_vals.get_cols());
    this->zero_val = _zero_val;

    size_t n_vals = this->rows*this->cols;
    this->vals.resize((n_vals));

    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    if (!nearly_equal<T>(_vals(i, j), _zero_val)) {
		this->vals.set(i*this->get_cols() + j, _vals(i, j));
	    }
	}
    }
}

// Initialize from a 2D vector
template<typename T>
SparseIntegerTypeMatrix<T>::SparseIntegerTypeMatrix(const std::vector<std::vector<T>> &rhs, const T& _zero_val) {
    this->resize_rows(rhs.size());
    this->resize_cols(rhs.at(0).size());
    this->zero_val = _zero_val;

    size_t n_vals = this->rows*this->cols;
    this->vals.resize(n_vals);

    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    if (!nearly_equal<T>(rhs[i][j], _zero_val)) {
		this->vals.set(i*this->get_cols() + j, rhs[i][j]);
	    }
	}
    }
}

// Copy constructor from contiguous 2D vector
template <typename T>
SparseIntegerTypeMatrix<T>::SparseIntegerTypeMatrix(const std::vector<T> &rhs, const size_t _rows, const size_t _cols, const T& _zero_val) {
    this->resize_rows(rhs.size());
    this->resize_cols(rhs.at(0).size());
    this->zero_val = _zero_val;

    size_t n_vals = this->rows*this->cols;
    this->vals.resize(n_vals);

    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    if (!nearly_equal<T>(rhs[i][j], _zero_val)) {
		this->vals.set(i*this->get_cols() + j, rhs[i][j]);
	    }
	}
    }
}

// Access individual elements
template <typename T>
T& SparseIntegerTypeMatrix<T>::operator()(size_t row, size_t col) {
    size_t address = row*this->get_cols() + col;
    if (this->vals[address] != this->zero_val) {
	return this->vals[address];
    }
    return zero_val;
}

// Access individual elements
template <typename T>
const T& SparseIntegerTypeMatrix<T>::operator()(size_t row, size_t col) const {
    size_t address = row*this->get_cols() + col;
    if (this->vals[address] != this->zero_val) {
	return this->vals[address];
    }
    return zero_val;
}

// TODO implement sparse matrix operators
// see https://www.geeksforgeeks.org/operations-sparse-matrices/ for reference

// In-place matrix-matrix addition
template<typename T>
SparseIntegerTypeMatrix<T>& SparseIntegerTypeMatrix<T>::operator+=(const Matrix<T>& rhs) {
    // seamat::SparseIntegerTypeMatrix<T>::operator+=
    //
    // Add values from rhs to the calling matrix in-place.
    //
    //   Input:
    //     `rhs`: Matrix to add, must have the same dimensions as the caller.
    //
    throw std::runtime_error("Sparse matrix operators have not been implemented.");
    return *this;
}

// In-place matrix-matrix subtraction
template<typename T>
SparseIntegerTypeMatrix<T>& SparseIntegerTypeMatrix<T>::operator-=(const Matrix<T>& rhs) {
    // seamat::SparseIntegerTypeMatrix<T>::operator-=
    //
    // Subtract values of rhs from the calling matrix in-place.
    //
    //   Input:
    //     `rhs`: Matrix to subtract, must have the same dimensions as the caller.
    //
    throw std::runtime_error("Sparse matrix operators have not been implemented.");

    return *this;
}

// In-place right multiplication
template<typename T>
SparseIntegerTypeMatrix<T>& SparseIntegerTypeMatrix<T>::operator*=(const Matrix<T>& rhs) {
    // seamat::SparseIntegerTypeMatrix<T>::operator*=
    //
    // Matrix right-multiplication of the caller with rhs in-place.
    //
    //   Input:
    //     `rhs`: Matrix to right multiply with,
    //            must have the same number of rows as the caller has columns.
    //
    throw std::runtime_error("Sparse matrix operators have not been implemented.");
    return *this;
}

// In-place left multiplication
template<typename T>
SparseIntegerTypeMatrix<T>& SparseIntegerTypeMatrix<T>::operator%=(const Matrix<T>& lhs) {
    // seamat::SparseIntegerTypeMatrix<T>::operator%=
    //
    // Matrix left-multiplication of the caller with lhs in-place.
    //
    //   Input:
    //     `lhs`: Matrix to left multiply with,
    //            must have the same number of columns as the caller has rows.
    //
    throw std::runtime_error("Sparse matrix operators have not been implemented.");
    return *this;
}

// In-place matrix-scalar addition
template<typename T>
SparseIntegerTypeMatrix<T>& SparseIntegerTypeMatrix<T>::operator+=(const T& scalar) {
    // seamat::SparseIntegerTypeMatrix<T>::operator+=
    //
    // In-place addition of a scalar to caller.
    //
    //   Input:
    //     `scalar`: Scalar value to add to all caller values.
    //
    this->zero_val += scalar;
    bm::bvector<>::enumerator en = this->vals.first();
    bm::bvector<>::enumerator en_end = this->vals.last();

    while (en < en_end) {
	*en += scalar;
	++en;
    }

    return *this;
}

// In-place matrix-scalar subtraction
template<typename T>
SparseIntegerTypeMatrix<T>& SparseIntegerTypeMatrix<T>::operator-=(const T& scalar) {
    // seamat::SparseIntegerTypeMatrix<T>::operator-=
    //
    // In-place subtraction of a scalar from the caller.
    //
    //   Input:
    //     `scalar`: Scalar value to subtract from all caller values.
    //
    this->zero_val -= scalar;
    bm::bvector<>::enumerator en = this->vals.first();
    bm::bvector<>::enumerator en_end = this->vals.last();

    while (en < en_end) {
	*en -= scalar;
	++en;
    }

    return *this;
}

// In-place matrix-scalar multiplication
template<typename T>
SparseIntegerTypeMatrix<T>& SparseIntegerTypeMatrix<T>::operator*=(const T& scalar) {
    // seamat::SparseIntegerTypeMatrix<T>::operator*=
    //
    // In-place multiplication of the caller with a scalar.
    //
    //   Input:
    //     `scalar`: Scalar value to multiply all caller values with.
    //
    // Handle special case where the whole matrix is multiplied by zero and becomes sparse.
    if (nearly_equal<T>(scalar, (T)0)) {
	this->vals.clear(true); // Free stored values
	this->vals.resize((size_t)this->get_rows*this->get_cols);
	this->zero_val = (T)0;
    } else {
	this->zero_val *= scalar;
	bm::bvector<>::enumerator en = this->vals.first();
	bm::bvector<>::enumerator en_end = this->vals.last();

	while (en < en_end) {
	    *en *= scalar;
	    ++en;
	}
    }

    return *this;
}

// In-place matrix-scalar division
template<typename T>
SparseIntegerTypeMatrix<T>& SparseIntegerTypeMatrix<T>::operator/=(const T& scalar) {
    // seamat::SparseIntegerTypeMatrix<T>::operator/=
    //
    // In-place division of the caller with a scalar.
    //
    //   Input:
    //     `scalar`: Scalar value to divide all caller values with.
    //
    if (nearly_equal<T>(scalar, (T)0))
	throw std::runtime_error("Math error: attempt to divide by zero.");

    this->zero_val *= scalar;
    bm::bvector<>::enumerator en = this->vals.first();
    bm::bvector<>::enumerator en_end = this->vals.last();

    while (en < en_end) {
	*en /= scalar;
	++en;
    }
}
}

#endif

// seamat: templatized matrix library
// https://github.com/tmaklin/seamat
//
// Copyright (C) 2021 Tommi Mäklin (tommi@maklin.fi)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef SEAMAT_SPARSE_MATRIX_CPP
#define SEAMAT_SPARSE_MATRIX_CPP

#include <cmath>
#include <stdexcept>

#include "Matrix.hpp"
#include "openmp_config.hpp"
#include "math_util.hpp"

namespace seamat {
template <typename T> class SparseMatrix : public Matrix<T> {
private:
    // TODO generic sparse matrix (instead of zero)
    //
    // Sparse matrix implemented in the compressed row storage (CRS) format.
    // See link below for reference.
    // https://netlib.org/linalg/html_templates/node91.html#SECTION00931100000000000000
    //
    std::vector<T> vals;
    std::vector<size_t> row_ptr;
    std::vector<size_t> col_ind;
    T zero_val;

    T* get_address(size_t row, size_t col);
    const T* get_address(size_t row, size_t col) const;

public:
    SparseMatrix() = default;
    ~SparseMatrix() = default;
    // Parameter constructor
    SparseMatrix(size_t _rows, size_t _cols, const T& _initial);
    // Initialize from a DenseMatrix
    SparseMatrix(const Matrix<T> &_vals, const T& _zero_val);
    // Initialize from a 2D vector
    SparseMatrix(const std::vector<std::vector<T>> &rhs, const T& _zero_val);
    // Copy constructor from contiguous 2D vector
    SparseMatrix(const std::vector<T> &rhs, const size_t _rows, const size_t _cols, const T& _zero_val);

    // Access individual elements
    T& operator()(size_t row, size_t col) override;
    const T& operator()(size_t row, size_t col) const override;

    // Mathematical operators
    // Matrix-matrix in-place summation and subtraction
    SparseMatrix<T>& operator+=(const Matrix<T>& rhs) override;
    SparseMatrix<T>& operator-=(const Matrix<T>& rhs) override;

    // In-place right multiplication
    SparseMatrix<T>& operator*=(const Matrix<T>& rhs) override;
    // In-place left multiplication
    SparseMatrix<T>& operator%=(const Matrix<T>& rhs) override;

    // Matrix-scalar, in-place
    SparseMatrix<T>& operator+=(const T& rhs) override;
    SparseMatrix<T>& operator-=(const T& rhs) override;
    SparseMatrix<T>& operator*=(const T& rhs) override;
    SparseMatrix<T>& operator/=(const T& rhs) override;
};

template<typename T>
T* SparseMatrix<T>::get_address(size_t row, size_t col) {
    // Returns the position of element (i, j) in this->vals
    size_t row_start = this->row_ptr[row];
    size_t row_end = this->row_ptr[row + 1];
    size_t nnz_row = row_end - row_start; // Number of non-zero elements in row `row`.
    if (nnz_row > 0) {
	for (size_t i = 0; i < nnz_row; ++i) {
	    size_t index = row_start + i;
	    if (this->col_ind[index] == col) {
		return &this->vals[index];
	    }
	}
    }
    return NULL;
}

template<typename T>
const T* SparseMatrix<T>::get_address(size_t row, size_t col) const {
    // Returns the position of element (i, j) in this->vals
    size_t row_start = this->row_ptr[row];
    size_t row_end = this->row_ptr[row + 1];
    size_t nnz_row = row_end - row_start; // Number of non-zero elements in row `row`.
    if (nnz_row > 0) {
	for (size_t i = 0; i < nnz_row; ++i) {
	    size_t index = row_start + i;
	    if (this->col_ind[index] == col) {
		return &this->vals[index];
	    }
	}
    }
    return NULL;
}

// Parameter constructor
template<typename T>
SparseMatrix<T>::SparseMatrix(size_t _rows, size_t _cols, const T& _initial) {
    // Initializes a dense SparseMatrix; use
    // SparseMatrix<T>::remove_nonzeros to sparsify the matrix after
    // filling it.
    this->resize_rows(_rows);
    this->resize_cols(_cols);

    uint64_t n_elements = _rows*_cols;
    this->vals = std::vector<T>(n_elements, _initial);
    this->col_ind.resize(n_elements, 0);
    this->row_ptr.resize(_rows, 0);
    for (size_t i = 0; i < _rows; ++i) {
	this->row_ptr[i + 1] = this->row_ptr[i] + _cols;
	for (size_t j = 0; j < _cols; ++j) {
	    this->col_ind[i*_cols + j] = j;
	}
    }
}
    
// Initialize from a DenseMatrix
template<typename T>
SparseMatrix<T>::SparseMatrix(const Matrix<T> &_vals, const T& _zero_val) {
    this->resize_rows(_vals.get_rows());
    this->resize_cols(_vals.get_cols());
    this->zero_val = _zero_val;

    uint64_t n_nonzero_elem = 0;
    for (size_t i = 0; i < _vals.get_rows(); ++i) {
	for (size_t j = 0; j < _vals.get_cols(); ++j) {
	    if (!nearly_equal<T>(_vals(i, j), this->zero_val)) {
		++n_nonzero_elem;
	    }
	}
    }

    this->row_ptr.resize(_vals.get_rows() + 1, 0);
    this->col_ind.resize(n_nonzero_elem, 0);
    this->vals.resize(n_nonzero_elem, _zero_val);

    size_t index = 0;
    for (size_t i = 0; i < this->get_rows(); ++i) {
	this->row_ptr[i + 1] = this->row_ptr[i];
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    const T &rhs_val = _vals(i, j);
	    if (!nearly_equal<T>(rhs_val, this->zero_val)) {
		++this->row_ptr[i + 1];
		this->col_ind[index] = j;
		this->vals[index] = rhs_val;
		++index;
	    }
	}
    }
}

// Initialize from a 2D vector
template<typename T>
SparseMatrix<T>::SparseMatrix(const std::vector<std::vector<T>> &rhs, const T& _zero_val) {
    this->resize_rows(rhs.size());
    this->resize_cols(rhs.at(0).size());
    this->zero_val = _zero_val;

    uint64_t n_nonzero_elem = 0;
    for (size_t i = 0; i < this->get_rows(); ++i) {
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    if (rhs[i][j] != this->zero_val) { // todo: floating point comparisons
		++n_nonzero_elem;
	    }
	}
    }

    this->row_ptr.resize(this->get_rows() + 1, 0);
    this->col_ind.resize(n_nonzero_elem, 0);
    this->vals.resize(n_nonzero_elem, _zero_val);

    size_t index = 0;
    for (size_t i = 0; i < this->get_rows(); ++i) {
	this->row_ptr[i + 1] = this->row_ptr[i];
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    if (rhs[i][j] != this->zero_val) { // todo: use std::nextafter for floating point comparisons
		++this->row_ptr[i + 1];
		this->col_ind[index] = j;
		this->vals[index] = rhs[i][j];
		++index;
	    }
	}
    }
}

// Copy constructor from contiguous 2D vector
template <typename T>
SparseMatrix<T>::SparseMatrix(const std::vector<T> &rhs, const size_t _rows, const size_t _cols, const T& _zero_val) {
    this->resize_rows(_rows);
    this->resize_cols(_cols);
    this->zero_val = _zero_val;

    uint64_t n_nonzero_elem = 0;
    for (size_t i = 0; i < _rows; ++i) {
	for (size_t j = 0; j < _cols; ++j) {
	    if (!nearly_equal<T>(rhs[i*_cols + j], this->zero_val)) {
		++n_nonzero_elem;
	    }
	}
    }

    this->row_ptr.resize(_rows + 1, 0);
    this->col_ind.resize(n_nonzero_elem, 0);
    this->vals.resize(n_nonzero_elem, _zero_val);

    size_t index = 0;
    for (size_t i = 0; i < this->get_rows(); ++i) {
	this->row_ptr[i + 1] = this->row_ptr[i];
	for (size_t j = 0; j < this->get_cols(); ++j) {
	    const T &rhs_val = rhs[i*_cols + j];
	    if (!nearly_equal<T>(rhs_val, this->zero_val)) {
		++this->row_ptr[i + 1];
		this->col_ind[index] = j;
		this->vals[index] = rhs_val;
		++index;
	    }
	}
    }
}

// Access individual elements
template <typename T>
T& SparseMatrix<T>::operator()(size_t row, size_t col) {
    T* address = this->get_address(row, col);
    if (address == NULL) {
	return this->zero_val;
    }
    return *address;
}

// Access individual elements (const)
template <typename T>
const T& SparseMatrix<T>::operator()(size_t row, size_t col) const {
    const T* address = this->get_address(row, col);
    if (address == NULL) {
	return this->zero_val;
    }
    return *address;
}


// TODO implement sparse matrix operators
// see https://www.geeksforgeeks.org/operations-sparse-matrices/ for reference

// In-place matrix-matrix addition
template<typename T>
SparseMatrix<T>& SparseMatrix<T>::operator+=(const Matrix<T>& rhs) {
    // seamat::SparseMatrix<T>::operator+=
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
SparseMatrix<T>& SparseMatrix<T>::operator-=(const Matrix<T>& rhs) {
    // seamat::SparseMatrix<T>::operator-=
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
SparseMatrix<T>& SparseMatrix<T>::operator*=(const Matrix<T>& rhs) {
    // seamat::SparseMatrix<T>::operator*=
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
SparseMatrix<T>& SparseMatrix<T>::operator%=(const Matrix<T>& lhs) {
    // seamat::SparseMatrix<T>::operator%=
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
SparseMatrix<T>& SparseMatrix<T>::operator+=(const T& scalar) {
    // seamat::SparseMatrix<T>::operator+=
    //
    // In-place addition of a scalar to caller.
    //
    //   Input:
    //     `scalar`: Scalar value to add to all caller values.
    //
    this->zero_val += scalar;
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->vals.size(); ++i) {
	this->vals[i] += scalar;
    }

    return *this;
}

// In-place matrix-scalar subtraction
template<typename T>
SparseMatrix<T>& SparseMatrix<T>::operator-=(const T& scalar) {
    // seamat::SparseMatrix<T>::operator-=
    //
    // In-place subtraction of a scalar from the caller.
    //
    //   Input:
    //     `scalar`: Scalar value to subtract from all caller values.
    //
    this->zero_val -= scalar;
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->vals.size(); ++i) {
	this->vals[i] -= scalar;
    }

    return *this;
}

// In-place matrix-scalar multiplication
template<typename T>
SparseMatrix<T>& SparseMatrix<T>::operator*=(const T& scalar) {
    // seamat::SparseMatrix<T>::operator*=
    //
    // In-place multiplication of the caller with a scalar.
    //
    //   Input:
    //     `scalar`: Scalar value to multiply all caller values with.
    //
    // Handle special case where the whole matrix is multiplied by zero and becomes sparse.
    if (nearly_equal<T>(scalar, (T)0)) {
	this->vals.clear();
	this->vals.shrink_to_fit();
	this->row_ptr.clear();
	this->row_ptr.shrink_to_fit();
	this->col_ind.clear();
	this->col_ind.shrink_to_fit();
	this->zero_val = (T)0;
    } else {
	this->zero_val *= scalar;
#pragma omp parallel for schedule(static)
	for (size_t i = 0; i < this->vals.size(); ++i) {
	    this->vals[i] *= scalar;
	}
    }

    return *this;
}

// In-place matrix-scalar division
template<typename T>
SparseMatrix<T>& SparseMatrix<T>::operator/=(const T& scalar) {
    // seamat::SparseMatrix<T>::operator/=
    //
    // In-place division of the caller with a scalar.
    //
    //   Input:
    //     `scalar`: Scalar value to divide all caller values with.
    //
    if (nearly_equal<T>(scalar, (T)0))
	throw std::runtime_error("Math error: attempt to divide by zero.");

    this->zero_val /= scalar;
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < this->vals.size(); ++i) {
	this->vals[i] /= scalar;
    }

    return *this;
}
}

#endif

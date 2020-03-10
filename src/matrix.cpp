#ifndef MATRIX_CPP
#define MATRIX_CPP
#include "matrix.hpp"

#include <cmath>

#if defined(MSWEEP_OPENMP_SUPPORT) && (MSWEEP_OPENMP_SUPPORT) == 1
#include <omp.h>
#include <algorithm>
#pragma omp declare reduction(vec_double_plus : std::vector<double> :	\
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
#endif

// Parameter Constructor
template<typename T>
Matrix<T>::Matrix(unsigned _rows, unsigned _cols, const T& _initial) {
  mat.resize(_rows);
#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < mat.size(); i++) {
    mat[i].resize(_cols, _initial);
  }
  rows = _rows;
  cols = _cols;
}

// Copy constructor
template<typename T>
Matrix<T>::Matrix(const Matrix<T>& rhs) {
  mat = rhs.mat;
  rows = rhs.get_rows();
  cols = rhs.get_cols();
}

// (Virtual) Destructor
template<typename T>
Matrix<T>::~Matrix() {}

// Resize a matrix
template<typename T>
void Matrix<T>::resize(const uint32_t new_rows, const uint32_t new_cols, const T initial) {
  if (new_rows != rows) {
    mat.resize(new_rows, std::vector<T>(new_cols, initial));
    rows = new_rows;
    cols = new_cols;
  }
  if (new_cols != cols) {
#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < new_rows; ++i) {
      mat[i].resize(new_cols, initial);
    }
    cols = new_cols;
  }
}


// Assignment Operator
template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& rhs) {
  if (&rhs == this)
    return *this;

  unsigned new_rows = rhs.get_rows();
  unsigned new_cols = rhs.get_cols();
  if (new_rows != rows || new_cols != cols) {
    resize(new_rows, new_cols, (T)0);
  }
#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < new_rows; i++) {
    for (unsigned j = 0; j < new_cols; j++) {
      mat[i][j] = rhs(i, j);
    }
  }
  return *this;
}

// Matrix-matrix addition
template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& rhs) const {
  Matrix result(this->rows, this->cols, 0.0);

#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < this->rows; i++) {
    for (unsigned j = 0; j < this->cols; j++) {
      result(i, j) = this->mat[i][j] + rhs(i,j);
    }
  }

  return result;
}

// In-place matrix-matrix addition
template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& rhs) {
#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < this->rows; i++) {
    for (unsigned j = 0; j < this->cols; j++) {
      this->mat[i][j] += rhs(i, j);
    }
  }

  return *this;
}

// Fill matrix with sum of two matrices
template <typename T>
void Matrix<T>::sum_fill(const Matrix<T>& rhs1, const Matrix<T>& rhs2) {
#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < this->rows; ++i) {
    for (unsigned j = 0; j < this->cols; ++j) {
      this->mat[i][j] = rhs1(i, j) + rhs2(i, j);
    }
  }
}

// Matrix-matrix subtraction
template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& rhs) const {
  Matrix result(this->rows, this->cols, 0.0);

#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < this->rows; i++) {
    for (unsigned j = 0; j < this->cols; j++) {
      result(i, j) = this->mat[i][j] - rhs(i, j);
    }
  }

  return result;
}

// In-place matrix-matrix subtraction
template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& rhs) {
#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < this->rows; i++) {
    for (unsigned j = 0; j < this->cols; j++) {
      this->mat[i][j] -= rhs(i, j);
    }
  }

  return *this;
}

// Matrix-matrix left multiplication
template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& rhs) const {
  Matrix result(this->rows, this->cols, 0.0);

#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < this->rows; i++) {
    for (unsigned j = 0; j < this->cols; j++) {
      for (unsigned k = 0; k < this->rows; k++) {
        result(i, j) += this->mat[i][k] * rhs(k, j);
      }
    }
  }

  return result;
}

// In-place matrix-matrix left multiplication
template<typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& rhs) {
  Matrix result = (*this) * rhs;
  (*this) = result;
  return *this;
}

// Transpose matrix
template<typename T>
Matrix<T> Matrix<T>::transpose() const {
  Matrix result(this->rows, this->cols, 0.0);

#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < this->rows; i++) {
    for (unsigned j = 0; j < this->cols; j++) {
      result(i, j) = this->mat[j][i];
    }
  }

  return result;
}

// In-place matrix-scalar addition
template<typename T>
Matrix<T>& Matrix<T>::operator+=(const T& rhs) {
#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < this->rows; i++) {
    for (unsigned j = 0; j < this->cols; j++) {
      this->mat[i][j] += rhs;
    }
  }

  return *this;
}

// In-place matrix-scalar subtraction
template<typename T>
Matrix<T>& Matrix<T>::operator-=(const T& rhs) {
#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < this->rows; i++) {
    for (unsigned j=0; j < this->cols; j++) {
      this->mat[i][j] -= rhs;
    }
  }

  return *this;
}

// In-place matrix-scalar multiplication
template<typename T>
Matrix<T>& Matrix<T>::operator*=(const T& rhs) {
#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < this->rows; ++i) {
    for (unsigned j = 0; j < this->cols; ++j) {
      this->mat[i][j] *= rhs;
    }
  }

  return *this;
}

// In-place matrix-scalar division
template<typename T>
Matrix<T>& Matrix<T>::operator/=(const T& rhs) {
#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < this->rows; ++i) {
    for (unsigned j = 0; j < this->cols; ++j) {
      this->mat[i][j] /= rhs;
    }
  }

  return *this;
}

// Matrix-vector right multiplication
template<typename T>
std::vector<T> Matrix<T>::operator*(const std::vector<T>& rhs) const {
  std::vector<T> result(rhs.size(), 0.0);

#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      result[i] += this->mat[i][j] * rhs[j];
    }
  }

  return result;
}

// Matrix-vector right multiplication, store result in arg
template<typename T>
void Matrix<T>::right_multiply(const std::vector<long unsigned>& rhs, std::vector<T>& result) const {
#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < this->rows; i++) {
    result[i] = 0.0;
    for (unsigned j = 0; j < this->cols; j++) {
      result[i] += this->mat[i][j] * rhs[j];
    }
  }
}

// log-space Matrix-vector right multiplication, store result in arg
template<typename T>
void Matrix<T>::exp_right_multiply(const std::vector<T>& rhs, std::vector<T>& result) const {
  std::fill(result.begin(), result.end(), 0.0);
#pragma omp parallel for schedule(static) reduction(vec_double_plus:result)
  for (unsigned i = 0; i < this->rows; i++) {
    for (unsigned j = 0; j < this->cols; j++) {
      result[i] += std::exp(this->mat[i][j] + rhs[j]);
    }
  }
}

// Specialized matrix-vector right multiplication
template<typename T>
std::vector<double> Matrix<T>::operator*(const std::vector<long unsigned>& rhs) const {
  std::vector<double> result(this->rows, 0.0);
  
#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < this->rows; i++) {
    for (unsigned j = 0; j < this->cols; j++) {
      result[i] += this->mat[i][j] * rhs[j];
    }
  }

  return result;
}

// Access individual elements
template<typename T>
T& Matrix<T>::operator()(unsigned row, unsigned col) {
  return this->mat[row][col];
}

// Access individual elements (const)
template<typename T>
const T& Matrix<T>::operator()(unsigned row, unsigned col) const {
  return this->mat[row][col];
}

// Access rows
template<typename T>
const std::vector<T>& Matrix<T>::get_row(unsigned row_id) const {
  return this->mat[row_id];
} 

// Access cols
template<typename T>
const std::vector<T>& Matrix<T>::get_col(unsigned col_id) const {
  std::vector<T> col(this->rows);
#pragma omp parallel for schedule(static)
  for (unsigned i = 0; i < this->rows; ++i) {
    col[i] = this->mat[i][col_id];
  }
  return col;
} 

// LogSumExp a Matrix column
template <typename T>
T Matrix<T>::log_sum_exp_col(unsigned col_id) const {
  // Note: this function accesses the elements rather inefficiently so
  // it shouldn't be parallellised here. However, the caller can
  // parallellize logsumexping multiple cols.
  T max_elem = 0;
  T sum = 0;
  for (unsigned i = 0; i < this->rows; ++i) {
    max_elem = (this->mat[i][col_id] > max_elem ? this->mat[i][col_id] : max_elem);
  }

  for (unsigned i = 0; i < this->rows; ++i) {
    sum += std::exp(this->mat[i][col_id] - max_elem);
  }
  return max_elem + std::log(sum);
}

// Get the number of rows of the matrix
template<typename T>
unsigned Matrix<T>::get_rows() const {
  return this->rows;
}

// Get the number of columns of the matrix
template<typename T>
unsigned Matrix<T>::get_cols() const {
  return this->cols;
}

#endif

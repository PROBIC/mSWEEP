#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

// Basic matrix structure and operations
// Implementation was done following the instructions at
// https://www.quantstart.com/articles/Matrix-Classes-in-C-The-Header-File
//
// **None of the operations validate the matrix sizes**

template <typename T> class Matrix {
 private:
  std::vector<std::vector<T> > mat;
  unsigned rows;
  unsigned cols;

 public:
  Matrix() = default;
  Matrix(unsigned _rows, unsigned _cols, const T& _initial);
  Matrix(const Matrix<T>& rhs);
  virtual ~Matrix();

  // Resize a matrix
  void resize(const uint32_t new_rows, const uint32_t new_cols, const T initial);

  // Operator overloading
  Matrix<T>& operator=(const Matrix<T>& rhs);

  // Mathematical operators
  // Matrix-matrix
  Matrix<T> operator+(const Matrix<T>& rhs) const;
  Matrix<T>& operator+=(const Matrix<T>& rhs);
  Matrix<T> operator-(const Matrix<T>& rhs) const;
  Matrix<T>& operator-=(const Matrix<T>& rhs);
  Matrix<T> operator*(const Matrix<T>& rhs) const;
  Matrix<T>& operator*=(const Matrix<T>& rhs);

  // Matrix-scalar, only in-place
  Matrix<T>& operator+=(const T& rhs);
  Matrix<T>& operator-=(const T& rhs);
  Matrix<T>& operator*=(const T& rhs);
  Matrix<T>& operator/=(const T& rhs);

  // Matrix-vector
  std::vector<T> operator*(const std::vector<T>& rhs) const;
  std::vector<double> operator*(const std::vector<long unsigned>& rhs) const;

  // Matrix-vector right multiplication, store result in arg
  void right_multiply(const std::vector<long unsigned>& rhs, std::vector<T>& result) const;
  void exp_right_multiply(const std::vector<T>& rhs, std::vector<T>& result) const;

  // Access elements
  T& operator()(unsigned row, unsigned col);
  const T& operator()(unsigned row, unsigned col) const;

  // Access rows and columns
  const std::vector<T>& get_row(unsigned row_id) const;
  const std::vector<T>& get_col(unsigned col_id) const;

  // LogSumExp a Matrix column
  T log_sum_exp_col(unsigned col_id) const;

  // Fill a matrix with the sum of two matrices
  void sum_fill(const Matrix<T>& rhs1, const Matrix<T>& rhs2);

  // Transpose
  Matrix<T> transpose() const;

  // Get number of rows or columns
  unsigned get_rows() const;
  unsigned get_cols() const;
};

#include "../src/matrix.cpp"

#endif

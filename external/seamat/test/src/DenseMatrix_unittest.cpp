// seamat: templatized matrix library
// https://github.com/tmaklin/seamat
//
// Copyright (C) 2021 Tommi Mäklin (tommi@maklin.fi)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#include "DenseMatrix_unittest.hpp"

TEST_F(DenseMatrixTest, ParameterConstructorWorks) {
    seamat::DenseMatrix<double> got(this->n_rows, this->n_cols, this->initial_val);
    EXPECT_EQ(this->n_rows, got.get_rows());
    EXPECT_EQ(this->n_cols, got.get_cols());
    for (uint32_t i = 0; i < this->n_rows; ++i) {
	for (uint32_t j = 0; j < this->n_cols; ++j) {
	    EXPECT_EQ(this->expected_mat[i*this->n_cols + j], got(i, j));
	}
    }
}

TEST_F(DenseMatrixTest, CopyConstructorWorks) {
    seamat::DenseMatrix<double> input(this->n_rows, this->n_cols, this->initial_val);
    seamat::DenseMatrix<double> got(input);
    EXPECT_EQ(this->n_rows, got.get_rows());
    EXPECT_EQ(this->n_cols, got.get_cols());
    for (uint32_t i = 0; i < this->n_rows; ++i) {
	for (uint32_t j = 0; j < this->n_cols; ++j) {
	    EXPECT_EQ(this->expected_mat[i*this->n_cols + j], got(i, j));
	}
    }
}

TEST_F(DenseMatrixTest, CopyConstructorFromVectorWorks) {
    seamat::DenseMatrix<double> got(this->expected_mat, this->n_rows, this->n_cols);
    EXPECT_EQ(this->n_rows, got.get_rows());
    EXPECT_EQ(this->n_cols, got.get_cols());
    for (uint32_t i = 0; i < this->n_rows; ++i) {
	for (uint32_t j = 0; j < this->n_cols; ++j) {
	    EXPECT_EQ(this->expected_mat[i*this->n_cols + j], got(i, j));
	}
    }
}

TEST_F(DenseMatrixTest, CopyConstructorFrom2DVectorWorks) {
    std::vector<std::vector<double>> input(3, std::vector<double>(7, 2.7));
    seamat::DenseMatrix<double> got(input);
    EXPECT_EQ(this->n_rows, got.get_rows());
    EXPECT_EQ(this->n_cols, got.get_cols());
    for (uint32_t i = 0; i < this->n_rows; ++i) {
	for (uint32_t j = 0; j < this->n_cols; ++j) {
	    EXPECT_EQ(this->expected_mat[i*this->n_cols + j], got(i, j));
	}
    }
}

TEST_F(DenseMatrixTest, ElementAssignmentWorks) {
    seamat::DenseMatrix<double> got(this->expected_mat, this->n_rows, this->n_cols);
    got(2, 5) = this->new_val;
    EXPECT_EQ(this->new_val, got(2, 5));
}

TEST_F(DenseMatrixTest, AssignmentOperatorWorks) {
    seamat::DenseMatrix<double> got(2, 5, 0.0);
    seamat::DenseMatrix<double> to_assign(this->expected_mat, this->n_rows, this->n_cols);
    got = to_assign;
    EXPECT_EQ(this->n_rows, got.get_rows());
    EXPECT_EQ(this->n_cols, got.get_cols());
    for (uint32_t i = 0; i < this->n_rows; ++i) {
	for (uint32_t j = 0; j < this->n_cols; ++j) {
	    EXPECT_EQ(this->expected_mat[i*this->n_cols + j], got(i, j));
	}
    }    
}

TEST_F(DenseMatrixTest, ResizeWorks) {
    seamat::DenseMatrix<double> got(2, 5, 0.0);
    got.resize(this->n_rows, this->n_cols, this->initial_val);
    EXPECT_EQ(this->n_rows, got.get_rows());
    EXPECT_EQ(this->n_cols, got.get_cols());
    for (uint32_t i = 0; i < 1; ++i) {
	for (uint32_t j = 0; j < 4; ++j) {
	    EXPECT_EQ(0.0, got(i, j));
	}
    }    
    for (uint32_t i = 1; i < this->n_rows; ++i) {
	for (uint32_t j = 4; j < this->n_cols; ++j) {
	    EXPECT_EQ(2.7, got(i, j));
	}
    }    
}

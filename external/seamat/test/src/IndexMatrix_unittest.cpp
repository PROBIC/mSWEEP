// seamat: templatized matrix library
// https://github.com/tmaklin/seamat
//
// Copyright (C) 2021 Tommi Mäklin (tommi@maklin.fi)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#include "IndexMatrix_unittest.hpp"

TEST_F(IndexMatrixTest, ValueConstructorFromDenseWorks) {
    std::vector<uint32_t> indices = { 1, 1, 1, 2, 0, 0, 0,
				      1, 1, 2, 2, 1, 0, 0,
				      2, 2, 0, 0, 0, 0, 0 };
    std::vector<double> vals = { 0.0, 1.0, 2.0,
				 0.0, 1.0, 2.0,
				 0.0, 1.0, 2.0 };
    seamat::IndexMatrix<double, uint32_t> got(seamat::DenseMatrix<double>(vals, 3, 3), seamat::DenseMatrix<uint32_t>(indices, 3, 7), false);
    EXPECT_EQ(this->n_rows, got.get_rows());
    EXPECT_EQ(this->n_cols, got.get_cols());
    for (uint32_t i = 0; i < this->n_rows; ++i) {
	for (uint32_t j = 0; j < this->n_cols; ++j) {
	    EXPECT_EQ(this->expected_mat[i*this->n_cols + j], got(i, j));
	}
    }
}

TEST_F(IndexMatrixTest, ValueConstructorFromSparseWorks) {
    std::vector<uint32_t> indices = { 1, 1, 1, 2, 0, 0, 0,
				      1, 1, 2, 2, 1, 0, 0,
				      2, 2, 0, 0, 0, 0, 0 };
    std::vector<double> vals = { 0.0, 1.0, 2.0,
				 0.0, 1.0, 2.0,
				 0.0, 1.0, 2.0 };
    seamat::IndexMatrix<double, uint32_t> got(seamat::DenseMatrix<double>(vals, 3, 3), seamat::DenseMatrix<uint32_t>(indices, 3, 7), true);
    EXPECT_EQ(this->n_rows, got.get_rows());
    EXPECT_EQ(this->n_cols, got.get_cols());
    for (uint32_t i = 0; i < this->n_rows; ++i) {
	for (uint32_t j = 0; j < this->n_cols; ++j) {
	    EXPECT_EQ(this->expected_mat[i*this->n_cols + j], got(i, j));
	}
    }
}

// seamat: templatized matrix library
// https://github.com/tmaklin/seamat
//
// Copyright (C) 2021 Tommi Mäklin (tommi@maklin.fi)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#ifndef SEAMAT_MATRIX_UNITTEST_HPP
#define SEAMAT_MATRIX_UNITTEST_HPP

#include "gtest/gtest.h"

#include "Matrix.hpp"

// Test parameter constructor
class MatrixOperatorsTest : public ::testing::Test {
    protected:
    void SetUp() override {
	this->n_rows = 3;
	this->n_cols = 7;

	// Set input matrix values
	this->lhs = { -1.250624, -0.2271989, -1.397814,   0.4060424, -0.1890402,  0.6845248, -0.102078,
		       1.273813,  0.3509934, -0.9456678, -0.859505,  -0.09857286, 0.3305294, -1.533049,
		      -0.5471792, 0.6436247, -0.3484785,  0.2681166, -0.05288153, 0.1891454, -0.4156959 };
	this->rhs = { -0.6097568, 1.146764,  -0.9064713, 0.8123275, 0.266325,  0.2133686,  1.348209,
		       1.16857,  -0.1501372, -0.3242086, 0.9417174, 0.3895409, 0.3988114, -1.589445,
		      -1.051972, -1.259515,  -0.9642974, 0.7157308, 2.730239,  0.3487662, -0.5616761 };
	this->expected = { -1.86038,   0.9195656, -2.304286, 1.21837,    0.07728473, 0.8978934,  1.24613,
			    2.442383,  0.2008561, -1.269876, 0.08221239, 0.290968,   0.7293408, -3.122493,
			   -1.599151, -0.6158905, -1.312776, 0.9838474,  2.677358,   0.5379115, -0.977372 };

	}
    void TearDown() override {
	this->n_rows = 0;
	this->n_cols = 0;
	this->lhs.clear();
	this->lhs.shrink_to_fit();
	this->rhs.clear();
	this->rhs.shrink_to_fit();
	this->expected.clear();
	this->expected.shrink_to_fit();
    }

    // Test inputs
    uint32_t n_rows;
    uint32_t n_cols;

    // Test input matrices
    std::vector<double> lhs;
    std::vector<double> rhs;

    // Expecteds
    std::vector<double> expected;
};

#endif

// mSWEEP: Estimate abundances of reference lineages in DNA sequencing reads.
//
// MIT License
//
// Copyright (c) 2023 Probabilistic Inference and Computational Biology group @ UH
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
#include "mSWEEP_io.hpp"

#include <sstream>
#include <exception>
#include <numeric>

#include "cxxio.hpp"

void WriteProbabilities(const seamat::DenseMatrix<double> &ec_probs, const std::vector<std::string> &cluster_indicators_to_string, std::ostream &of) {
  // Write the probability matrix to a file.
  if (of.good()) {
    of << "ec_id" << '\t';
    size_t n_rows = ec_probs.get_rows();
    size_t n_cols = ec_probs.get_cols();
    for (uint32_t i = 0; i < n_rows; ++i) {
      of << cluster_indicators_to_string[i];
      of << (i < n_rows - 1 ? '\t' : '\n');
    }
    for (uint32_t i = 0; i < n_cols; ++i) {
      of << i << '\t';
      for (uint32_t j = 0; j < n_rows; ++j) {
	of << std::exp(ec_probs(j, i));
	of << (j < n_rows - 1 ? '\t' : '\n');
      }
    }
    of << std::endl;
    of.flush();
  } else {
    throw std::runtime_error("Can't write to probs file.");
  }
}

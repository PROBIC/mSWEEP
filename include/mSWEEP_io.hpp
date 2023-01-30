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
#ifndef MSWEEP_MSWEEP_IO_HPP
#define MSWEEP_MSWEEP_IO_HPP

#include <string>
#include <vector>
#include <fstream>
#include <cstddef>

#include "telescope.hpp"
#include "Matrix.hpp"

#include "Reference.hpp"
#include "likelihood.hpp"

// Read functions
//// Read pseudoalignments
void ReadPseudoalignments(const std::vector<std::string> &alignment_paths, const std::string &themisto_merge_mode, const Reference &reference, std::unique_ptr<telescope::Alignment> &aln, LL_WOR21<double, uint16_t>* log_likelihoods);

// Write functions
//// Write likelihoods
void WriteLikelihood(const seamat::DenseMatrix<double> &ll_mat, const std::vector<double> &log_ec_counts, const uint32_t n_groups, std::ostream &of);
void WriteLikelihoodBitSeq(const seamat::DenseMatrix<double> &ll_mat, const std::vector<double> &log_ec_counts, const uint32_t n_groups, std::ostream &of);

//// Write read-reference probabilities
void WriteProbabilities(const seamat::DenseMatrix<double> &ec_probs, const std::vector<std::string> &cluster_indicators_to_string, std::ostream &of);

#endif

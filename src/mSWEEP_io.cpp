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

void ReadPseudoalignments(const std::vector<std::string> &alignment_paths, const std::string &themisto_merge_mode, const Reference &reference, std::unique_ptr<telescope::Alignment> &aln, LL_WOR21<double, uint16_t>* log_likelihoods) {
  size_t n_files = alignment_paths.size();
  std::vector<cxxio::In> infiles;
  infiles.reserve(n_files);
  std::vector<std::istream*> strands(n_files);
  if (n_files > 0) {
    for (size_t i = 0; i < n_files; ++i) {
      infiles.emplace_back(cxxio::In(alignment_paths[i]));
      strands[i] = &infiles[i].stream();
    }
  } else {
    strands.emplace_back(&std::cin);
  }
  telescope::read::ThemistoGrouped(telescope::get_mode(themisto_merge_mode), reference.get_n_refs(), reference.get_group_indicators(0), strands, aln);

  try {
    // Use the alignment data to populate the log_likelihoods matrix.
    log_likelihoods->from_grouped_alignment(*aln, reference.get_grouping(0));
  } catch (std::exception &e) {
    throw e;
  }
}

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

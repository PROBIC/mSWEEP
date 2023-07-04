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
#include "Sample.hpp"

#include <exception>

namespace mSWEEP {
void ConstructSample(const telescope::Alignment &alignment, const size_t bootstrap_iters, const size_t bootstrap_count, const size_t bootstrap_seed, const bool bin_reads, std::unique_ptr<Sample> &sample) {
  // Wrapper for determining which Sample type to construct.
  // Initialize Sample depending on how the alignment needs to be processed.
  if (bootstrap_iters > 0) {
    // Bootstrap mode
    bool count_provided = bootstrap_count > 0;
    if (count_provided && bin_reads) {
      sample.reset(new BinningBootstrap(alignment, bootstrap_iters, bootstrap_count, bootstrap_seed));
    } else if (count_provided && !bin_reads) {
      sample.reset(new BootstrapSample(alignment, bootstrap_iters, bootstrap_iters, bootstrap_seed));
    } else if (bin_reads) {
      sample.reset(new BinningBootstrap(alignment, bootstrap_iters, bootstrap_seed));
    } else {
      sample.reset(new BootstrapSample(alignment, bootstrap_iters, bootstrap_seed));
    }
  } else if (bin_reads) {
    sample.reset(new BinningSample(alignment));
  } else {
    sample.reset(new PlainSample(alignment));
  }
}

void Sample::count_alignments(const telescope::Alignment &alignment) {
  // Count the number of aligned reads and store in counts_total
  size_t aln_counts_total = 0;
#pragma omp parallel for schedule(static) reduction(+:aln_counts_total)
  for (size_t i = 0; i < alignment.n_ecs(); ++i) {
    aln_counts_total += alignment.reads_in_ec(i);
  }
  this->counts_total = aln_counts_total;
  this->n_reads = alignment.n_reads();
}

void Sample::write_probs(const std::vector<std::string> &cluster_indicators_to_string, std::ostream *of) {
  // Write the probability matrix to a file.
  if (of->good()) {
    *of << "ec_id" << '\t';
    size_t n_rows = this->ec_probabilities.get_rows();
    size_t n_cols = this->ec_probabilities.get_cols();
    for (size_t i = 0; i < n_rows; ++i) {
      *of << cluster_indicators_to_string[i];
      *of << (i < n_rows - 1 ? '\t' : '\n');
    }
    for (size_t i = 0; i < n_cols; ++i) {
	*of << i << '\t';
	for (size_t j = 0; j < n_rows; ++j) {
	  *of << std::exp(this->ec_probabilities(j, i));
	  *of << (j < n_rows - 1 ? '\t' : '\n');
	}
    }
    *of << std::endl;
      of->flush();
  } else {
    throw std::runtime_error("Can't write to probs file.");
  }
}

}

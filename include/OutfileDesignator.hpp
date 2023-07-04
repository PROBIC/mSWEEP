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
// mSWEEP::OutfileDesignator
// This file implements the OutfileDesignator class which handles setting the filenames
// correctly and providing the output streams.
//
#ifndef MSWEEP_OUTFILEDESIGNATOR_HPP
#define MSWEEP_OUTFILEDESIGNATOR_HPP

#include <string>
#include <cstddef>
#include <fstream>

#include "bxzstr.hpp"
#include "cxxio.hpp"

namespace mSWEEP {
class OutfileDesignator {
private:
  // Was printing abundances requested?
  bool printing;

  // Compression info for writing probs/likelihoods/bins
  bool compress;
  bxz::Compression type;
  int compression_level;
  std::string extension;

  // Filename info
  std::string prefix;
  size_t n_groupings;
  size_t current_grouping;

  // Out stream is wrapped in cxxio for checking writability etc.
  cxxio::Out of;

  // Helper for opening compressed or plain output.
  void open(std::string &filename);

public:
  // Throwing constructor
  OutfileDesignator(std::string _prefix, size_t _n_groupings, std::string _compress, int _level);

  // Likelihoods output file, format should be "bitseq" or anything else
  std::ostream* likelihoods(const std::string &format);

  // mGEMS bins, note these ignore the grouping (groups should have different names).
  std::ostream* bin(const std::string &name);

  // Probs
  std::ostream* probs();

  // Abundances, ignores compression
  std::ostream* abundances();

  // Increment the grouping counter
  void next_grouping();

};

}

#endif

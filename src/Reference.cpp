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
#include "Reference.hpp"

std::unique_ptr<Reference> ConstructAdaptiveReference(std::istream *in, const char delimiter) {
  std::vector<std::string> group_indicators;
  std::set<std::string> group_names;
  if (in->good()) {
    std::string line;
    while (std::getline(*in, line)) {
      group_indicators.emplace_back(line);
      group_names.insert(line);
    }
  } else {
    throw std::runtime_error("Could not read cluster indicators.");
  }
  size_t n_groups = group_names.size();

  std::unique_ptr<Reference> ret;
  if (n_groups <= std::numeric_limits<uint8_t>::max()) {
    ret.reset(new AdaptiveReference<uint8_t>(group_indicators, delimiter));
  } else if (n_groups <= std::numeric_limits<uint16_t>::max()) {
    ret.reset(new AdaptiveReference<uint16_t>(group_indicators, delimiter));
  } else if (n_groups <= std::numeric_limits<uint32_t>::max()) {
    ret.reset(new AdaptiveReference<uint32_t>(group_indicators, delimiter));
  } else {
    ret.reset(new AdaptiveReference<uint64_t>(group_indicators, delimiter));
  }
  return ret;
}

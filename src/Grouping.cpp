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
#include "Grouping.hpp"

std::unique_ptr<Grouping> ConstructAdaptive(const std::vector<std::string> &indicators) {
  std::unordered_map<std::string, size_t> name_to_size;
  for (size_t i = 0; i < indicators.size(); ++i) {
    if (name_to_size.find(indicators[i]) == name_to_size.end()) {
      name_to_size[indicators[i]] = 0;
    }
    ++name_to_size[indicators[i]];
  }

  size_t max_size = 0;
  for (auto kv : name_to_size) {
    max_size = (kv.second > max_size ? kv.second : max_size);
  }

  size_t n_groups = name_to_size.size();

  std::unique_ptr<Grouping> ret;
  if (max_size <= std::numeric_limits<uint8_t>::max()) {
    if (n_groups <= std::numeric_limits<uint8_t>::max()) {
      ret.reset(new AdaptiveGrouping<uint8_t, uint8_t>(indicators));
    } else if (n_groups <= std::numeric_limits<uint16_t>::max()) {
      ret.reset(new AdaptiveGrouping<uint8_t, uint16_t>(indicators));
    } else if (n_groups <= std::numeric_limits<uint32_t>::max()) {
      ret.reset(new AdaptiveGrouping<uint8_t, uint32_t>(indicators));
    } else {
      ret.reset(new AdaptiveGrouping<uint8_t, uint64_t>(indicators));
    }
  } else if (max_size <= std::numeric_limits<uint16_t>::max()) {
    if (n_groups <= std::numeric_limits<uint8_t>::max()) {
      ret.reset(new AdaptiveGrouping<uint8_t, uint8_t>(indicators));
    } else if (n_groups <= std::numeric_limits<uint16_t>::max()) {
      ret.reset(new AdaptiveGrouping<uint8_t, uint16_t>(indicators));
    } else if (n_groups <= std::numeric_limits<uint32_t>::max()) {
      ret.reset(new AdaptiveGrouping<uint8_t, uint32_t>(indicators));
    } else {
      ret.reset(new AdaptiveGrouping<uint8_t, uint64_t>(indicators));
    }
  } else if (max_size <= std::numeric_limits<uint32_t>::max()) {
    if (n_groups <= std::numeric_limits<uint8_t>::max()) {
      ret.reset(new AdaptiveGrouping<uint8_t, uint8_t>(indicators));
    } else if (n_groups <= std::numeric_limits<uint16_t>::max()) {
      ret.reset(new AdaptiveGrouping<uint8_t, uint16_t>(indicators));
    } else if (n_groups <= std::numeric_limits<uint32_t>::max()) {
      ret.reset(new AdaptiveGrouping<uint8_t, uint32_t>(indicators));
    } else {
      ret.reset(new AdaptiveGrouping<uint8_t, uint64_t>(indicators));
    }
  } else {
    if (n_groups <= std::numeric_limits<uint8_t>::max()) {
      ret.reset(new AdaptiveGrouping<uint8_t, uint8_t>(indicators));
    } else if (n_groups <= std::numeric_limits<uint16_t>::max()) {
      ret.reset(new AdaptiveGrouping<uint8_t, uint16_t>(indicators));
    } else if (n_groups <= std::numeric_limits<uint32_t>::max()) {
      ret.reset(new AdaptiveGrouping<uint8_t, uint32_t>(indicators));
    } else {
      ret.reset(new AdaptiveGrouping<uint8_t, uint64_t>(indicators));
    }
  }
  return ret;
}

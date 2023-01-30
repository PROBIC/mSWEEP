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

#include <sstream>
#include <exception>

#include "tools/matchfasta.hpp"

void Reference::add_sequence(const std::string &indicator_s, const uint16_t grouping_id) {
  if (grouping_id == 0) {
    this->n_refs += 1;
  }
  this->groups_indicators[grouping_id].emplace_back(this->groupings[grouping_id]->get_id(indicator_s));
}

void Reference::read_from_file(std::istream &indicator_file, const char delimiter) {
  if (indicator_file.good()) {
    std::vector<std::vector<std::string>> group_indicators;
    std::string indicator_s;
    while (std::getline(indicator_file, indicator_s)) {
      std::stringstream indicators(indicator_s);
      std::string indicator;
      uint16_t grouping_id = 0;
      while (std::getline(indicators, indicator, delimiter)) {
	if (grouping_id >= this->n_groupings) {
	  this->groups_indicators.emplace_back(std::vector<uint32_t>());
	  group_indicators.emplace_back(std::vector<std::string>());
	  ++this->n_groupings;
	}
	group_indicators[grouping_id].emplace_back(indicator);
	++grouping_id;
      }
    }

    for (size_t i = 0; i < n_groupings; ++i) {
      this->groupings.emplace_back(ConstructAdaptive(group_indicators[i]));
      for (size_t j = 0; j < group_indicators[i].size(); ++j) {
	this->add_sequence(group_indicators[i][j], i);
      }
    }

  } else {
    throw std::runtime_error("Could not read cluster indicators.");
  }
  if (this->n_refs == 0) {
    throw std::runtime_error("The grouping contains 0 reference sequences");
  }
}

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
#ifndef MSWEEP_REFERENCE_HPP
#define MSWEEP_REFERENCE_HPP

#include <cstddef>
#include <vector>
#include <string>
#include <memory>
#include <sstream>
#include <exception>

#include "Grouping.hpp"

class Reference {
public:
  // Getters to access the groupings
  virtual size_t n_groups(const size_t grouping_id) const =0;
  virtual const std::vector<std::string>& group_names(const size_t grouping_id) const =0;

  // Getters
  virtual const Grouping& get_grouping(const size_t grouping_id) const =0;
  virtual size_t get_n_refs() const =0;
  virtual size_t get_n_groupings() const =0;
};

template <typename T>
class AdaptiveReference : public Reference {
private:
  size_t n_refs = 0;
  size_t n_groupings = 0;

  std::vector<std::vector<T>> groups_indicators;
  std::vector<std::unique_ptr<Grouping>> groupings;

  void add_sequence(const std::string &seq_name, const size_t grouping_id) {
    if (grouping_id == 0) {
      this->n_refs += 1;
    }
    this->groups_indicators[grouping_id].emplace_back(this->groupings[grouping_id]->get_id(seq_name));
  }

public:
  AdaptiveReference(const std::vector<std::string> &indicator_lines, const char delimiter = '\t') {
    std::vector<std::vector<std::string>> group_indicators;
    for (size_t i = 0; i < indicator_lines.size(); ++i) {
      std::string indicator_s;
      std::stringstream indicators(indicator_lines[i]);
      std::string indicator;
      size_t grouping_id = 0;
      while (std::getline(indicators, indicator, delimiter)) {
	if (grouping_id >= this->n_groupings) {
	  this->groups_indicators.emplace_back(std::vector<T>());
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
    if (this->n_refs == 0) {
      throw std::runtime_error("The grouping contains 0 reference sequences");
    }
  }

  // Getters to access the groupings
  size_t n_groups(const size_t grouping_id) const override { return (size_t)(*this->groupings[grouping_id]).get_n_groups(); }
  const std::vector<std::string>& group_names(const size_t grouping_id) const override { return (*this->groupings[grouping_id]).get_names(); };

  template <typename U>
  const std::vector<U>& group_sizes(const size_t grouping_id) const {
    return static_cast<const AdaptiveGrouping<U, T>*>(&(*this->groupings[grouping_id]))->get_sizes();
  }

  // Getters
  const Grouping& get_grouping(const size_t grouping_id) const override { return (*this->groupings[grouping_id]); };
  const std::vector<T>& get_group_indicators(const size_t grouping_id) const { return this->groups_indicators[grouping_id]; };
  size_t get_n_refs() const override { return (size_t)this->n_refs; };
  size_t get_n_groupings() const override { return (size_t)this->n_groupings; };

};

std::unique_ptr<Reference> ConstructAdaptiveReference(std::istream *in, const char delimiter = '\t');

#endif

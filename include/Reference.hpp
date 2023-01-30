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
#include <fstream>
#include <memory>

#include "Grouping.hpp"

class Reference {
private:
  uint32_t n_refs = 0;
  uint16_t n_groupings = 0;

  std::vector<std::vector<uint32_t>> groups_indicators;
  std::vector<std::unique_ptr<Grouping>> groupings;

  void add_sequence(const std::string &seq_name, const uint16_t grouping_id);

public:
  void read_from_file(std::istream &indicator_file, const char delimiter = '\t');
  void match_with_fasta(const char delimiter, std::istream &groups_file, std::istream &fasta_file);

  // Getters to access the groupings
  uint32_t n_groups(const size_t grouping_id) const { return (*this->groupings[grouping_id]).get_n_groups(); }
  const std::vector<std::string>& group_names(const size_t grouping_id) const { return (*this->groupings[grouping_id]).get_names(); };
  template <typename T, typename V>
  const std::vector<T>& group_sizes(const size_t grouping_id) const {
    return static_cast<const AdaptiveGrouping<T, V>*>(&(*this->groupings[grouping_id]))->get_sizes();
  }

  // Getters
  const Grouping& get_grouping(const uint16_t grouping_id) const { return (*this->groupings[grouping_id]); };
  const std::vector<uint32_t>& get_group_indicators(const uint16_t grouping_id) const { return this->groups_indicators[grouping_id]; };
  uint32_t get_n_refs() const { return this->n_refs; };
  uint16_t get_n_groupings() const { return this->n_groupings; };

};

#endif

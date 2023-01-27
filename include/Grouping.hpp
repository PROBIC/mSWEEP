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
#ifndef MSWEEP_GROUPING_HPP
#define MSWEEP_GROUPING_HPP

#include <cstddef>
#include <vector>
#include <string>
#include <unordered_map>
#include <array>

class Grouping {
private:
  uint32_t n_groups = 0;

  std::vector<uint16_t> sizes;
  std::vector<std::string> names;
  std::unordered_map<std::string, uint32_t> name_to_id;

  // Adds group_name to this Grouping
  void add_group(const std::string &group_name);

public:
  // Increase the count of sequences assigned to group_name, and add
  // group_name to this grouping if it does not exist yet.
  void add_sequence(const std::string &group_name);

  // Calculates the group-specific likelihood parameters that depend on the group sizes
  std::vector<std::array<double, 2>> bb_parameters(const double bb_constants[2]) const;

  // Find the numeric id of a group by its name
  uint32_t get_id(const std::string &name) const { return this->name_to_id.at(name); };

  // Find size of the largest group
  size_t max_group_size() const {
    size_t max_size = 0;
    for (size_t i = 0; i < this->n_groups; ++i)
      max_size = (this->sizes[i] > max_size ? this->sizes[i] : max_size);
    return max_size;
  }

  // Getters
  const std::vector<std::string>& get_names() const { return this->names; };
  const std::vector<uint16_t>& get_sizes() const { return this->sizes; };
  uint32_t get_n_groups() const { return this->n_groups; };
};

#endif

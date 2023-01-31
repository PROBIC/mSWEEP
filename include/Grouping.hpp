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

#include <string>
#include <cstddef>
#include <vector>
#include <unordered_map>
#include <memory>

class Grouping {
public:
  // Increase the count of sequences assigned to group_name, and add
  // group_name to this grouping if it does not exist yet.
  virtual void add_sequence(const std::string &group_name) =0;

  // Find the numeric id of a group by its name
  virtual size_t get_id(const std::string &name) const =0;

  // Find size of the largest group
  virtual size_t max_group_size() const =0;

  // Getters
  virtual const std::vector<std::string>& get_names() const =0;
  virtual size_t get_n_groups() const =0;
};

template <typename T, typename V>
class AdaptiveGrouping : public Grouping {
private:
  V n_groups = 0;

  std::vector<T> sizes;
  std::vector<std::string> names;
  std::unordered_map<std::string, V> name_to_id;

  // Adds group_name to this Grouping
  void add_group(const std::string &group_name) {
    this->name_to_id[group_name] = this->name_to_id.size(); // Newest always has id equal to size.
    this->names.emplace_back(group_name);
    this->sizes.emplace_back(0);
    this->n_groups += 1;
  }

public:
  AdaptiveGrouping(const std::vector<std::string> &indicators) {
    for (size_t i = 0; i < indicators.size(); ++i) {
      this->add_sequence(indicators[i]);
    }
  }

  // Increase the count of sequences assigned to group_name, and add
  // group_name to this grouping if it does not exist yet.
  void add_sequence(const std::string &group_name) override {
    if (this->name_to_id.find(group_name) == this->name_to_id.end()) {
      this->add_group(group_name);
    }
    this->sizes[this->name_to_id[group_name]] += 1;
  }

  // Find the numeric id of a group by its name
  size_t get_id(const std::string &name) const override { return (size_t)this->name_to_id.at(name); };

  // Find size of the largest group
  size_t max_group_size() const override {
    size_t max_size = 0;
    for (size_t i = 0; i < this->n_groups; ++i)
      max_size = (this->sizes[i] > max_size ? this->sizes[i] : max_size);
    return max_size;
  }

  // Getters
  const std::vector<std::string>& get_names() const override { return this->names; };
  const std::vector<T>& get_sizes() const { return this->sizes; };
  size_t get_n_groups() const override { return (size_t)this->n_groups; };
};

std::unique_ptr<Grouping> ConstructAdaptive(const std::vector<std::string> &indicators);

#endif

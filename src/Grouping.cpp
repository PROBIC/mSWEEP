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

std::vector<std::array<double, 2>> Grouping::bb_parameters(const double bb_constants[2]) const {
  std::vector<std::array<double, 2>> bb_params(this->n_groups);
  for (size_t i = 0; i < this->n_groups; ++i) {
    double e = this->sizes[i]*bb_constants[0];
    double phi = 1.0/(this->sizes[i] - e + bb_constants[1]);
    double beta = phi*(this->sizes[i] - e);
    double alpha = (e*beta)/(this->sizes[i] - e);
    bb_params[i] = std::array<double, 2>{ { alpha, beta } };
  }
  return bb_params;
}

void Grouping::add_group(const std::string &group_name) {
  this->name_to_id[group_name] = this->name_to_id.size(); // Newest always has id equal to size.
  this->names.emplace_back(group_name);
  this->sizes.emplace_back(0);
  this->n_groups += 1;
}

void Grouping::add_sequence(const std::string &group_name) {
  if (this->name_to_id.find(group_name) == this->name_to_id.end()) {
    this->add_group(group_name);
  }
  this->sizes[this->name_to_id[group_name]] += 1;
}

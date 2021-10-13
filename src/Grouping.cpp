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
  this->names.emplace_back(group_name);
  this->sizes.emplace_back(0);
  this->n_groups += 1;
}

void Grouping::add_sequence(const std::string &group_name) {
  this->sizes[this->name_to_id[group_name]] += 1;
}

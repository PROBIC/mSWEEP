#include "Grouping.hpp"

void Grouping::calculate_bb_parameters(double params[2]) {
  for (size_t i = 0; i < this->n_groups; ++i) {
    double e = this->sizes[i]*params[0];
    double phi = 1.0/(this->sizes[i] - e + params[1]);
    double beta = phi*(this->sizes[i] - e);
    double alpha = (e*beta)/(this->sizes[i] - e);
    this->bb_params.emplace_back(std::array<double, 2>{ { alpha, beta } });
  }
}
#include "Reference.hpp"

void Reference::calculate_bb_parameters(double params[2]) {
  for (size_t i = 0; i < this->grouping.n_groups; ++i) {
    double e = this->grouping.sizes[i]*params[0];
    double phi = 1.0/(this->grouping.sizes[i] - e + params[1]);
    double beta = phi*(this->grouping.sizes[i] - e);
    double alpha = (e*beta)/(this->grouping.sizes[i] - e);
    this->grouping.bb_params.emplace_back(std::array<double, 2>{ { alpha, beta } });
  }
}

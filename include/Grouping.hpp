#ifndef MSWEEP_GROUPING_HPP
#define MSWEEP_GROUPING_HPP

#include <vector>
#include <array>

class Grouping {
public:
  std::vector<uint32_t> indicators;
  std::vector<uint16_t> sizes;
  std::vector<std::array<double, 2>> bb_params;
  uint32_t n_groups;

  void calculate_bb_parameters(double params[2]);
};

#endif

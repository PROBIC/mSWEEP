#ifndef MSWEEP_GROUPING_HPP
#define MSWEEP_GROUPING_HPP

#include <vector>
#include <array>
#include <string>

class Grouping {
public:
  uint32_t n_groups;

  std::vector<uint32_t> indicators;
  std::vector<uint16_t> sizes;
  std::vector<std::string> names;

  std::vector<std::array<double, 2>> bb_params;

  void calculate_bb_parameters(double params[2]);
};

#endif

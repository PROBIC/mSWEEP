#ifndef REFERENCE_H
#define REFERENCE_H

#include <vector>
#include <array>

struct Grouping {
  std::vector<uint32_t> indicators;
  std::vector<uint16_t> sizes;
  std::vector<std::array<double, 2>> bb_params;
  uint32_t n_groups;
};

class Reference {
public:
  Grouping grouping;
  std::vector<std::string> group_names;
  uint32_t n_refs;
  
  void calculate_bb_parameters(double params[2]);
};

#endif

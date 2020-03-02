#ifndef REFERENCE_H
#define REFERENCE_H

#include <vector>
#include <array>

struct Grouping {
  std::vector<unsigned short> indicators;
  std::vector<short unsigned> sizes;
  std::vector<std::array<double, 2>> bb_params;
  unsigned short n_groups;
};

class Reference {
public:
  Grouping grouping;
  std::vector<std::string> group_names;
  unsigned n_refs;
  
  void calculate_bb_parameters(double params[2]);
};

#endif

#ifndef MSWEEP_GROUPING_HPP
#define MSWEEP_GROUPING_HPP

#include <vector>
#include <array>
#include <string>
#include <unordered_map>

class Grouping {
public:
  uint32_t n_groups = 0;

  std::vector<uint16_t> sizes;
  std::vector<std::string> names;
  std::unordered_map<std::string, uint32_t> name_to_id;

  void add_group(const std::string &group_name);
  void add_sequence(const std::string &group_name);

  // Calculates the group-specific likelihood parameters that depend on the group sizes
  std::vector<std::array<double, 2>> bb_parameters(const double bb_constants[2]) const;
};

#endif

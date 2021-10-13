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

  std::vector<std::array<double, 2>> bb_params;

  void calculate_bb_parameters(double params[2]);
  void add_group(const std::string &group_name);
  void add_sequence(const std::string &group_name);
};

#endif

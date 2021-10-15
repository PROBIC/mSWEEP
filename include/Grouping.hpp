#ifndef MSWEEP_GROUPING_HPP
#define MSWEEP_GROUPING_HPP

#include <vector>
#include <array>
#include <string>
#include <unordered_map>

class Grouping {
private:
  uint32_t n_groups = 0;

  std::vector<uint16_t> sizes;
  std::vector<std::string> names;
  std::unordered_map<std::string, uint32_t> name_to_id;

public:
  void add_group(const std::string &group_name);
  void add_sequence(const std::string &group_name);

  // Calculates the group-specific likelihood parameters that depend on the group sizes
  std::vector<std::array<double, 2>> bb_parameters(const double bb_constants[2]) const;

  // Find the numeric id of a group by its name
  uint32_t get_id(const std::string &name) const { return this->name_to_id.at(name); };

  // Getters
  const std::vector<std::string>& get_names() const { return this->names; };
  const std::vector<uint16_t>& get_sizes() const { return this->sizes; };
  uint32_t get_n_groups() const { return this->n_groups; };
};

#endif

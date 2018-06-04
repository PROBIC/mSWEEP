#ifndef READ_BITFIELD_H
#define READ_BITFIELD_H

#include <memory>
#include <string>
#include <unordered_map>

#include "matrix.hpp"
#include "Sample.hpp"
#include "Reference.hpp"

void ReadClusterIndicators(std::string &indicators_file, Reference &reference);
void ReadBitfield(std::vector<std::string> &kallisto_files, unsigned n_refs, std::vector<Sample> &batch);
void VerifyGrouping(std::string &run_info, unsigned n_refs);

#endif

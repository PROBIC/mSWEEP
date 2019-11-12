#ifndef READ_BITFIELD_H
#define READ_BITFIELD_H

#include <memory>
#include <string>
#include <unordered_map>
#include <fstream>

#include "matrix.hpp"
#include "Sample.hpp"
#include "Reference.hpp"
#include "zstr.hpp"
#include "KallistoFiles.hpp"

void ReadClusterIndicators(std::istream &indicators_file, Reference &reference);
void ReadBitfield(KallistoFiles &kallisto_files, unsigned n_refs, std::vector<Sample> &batch);
void VerifyGrouping(std::istream &run_info, unsigned n_refs);

#endif

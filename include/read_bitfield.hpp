#ifndef MSWEEP_READ_BITFIELD_HPP
#define MSWEEP_READ_BITFIELD_HPP

#include <vector>
#include <memory>
#include <string>

#include "matrix.hpp"
#include "Sample.hpp"
#include "Reference.hpp"
#include "KallistoFiles.hpp"

void ReadBitfield(KallistoFiles &kallisto_files, unsigned n_refs, std::vector<std::unique_ptr<Sample>> &batch, Reference &reference, bool bootstrap_mode);
void ReadBitfield(const std::string &tinfile1, const std::string &tinfile2, const std::string &themisto_mode, const bool bootstrap_mode, const unsigned n_refs, std::vector<std::unique_ptr<Sample>> &batch);

#endif

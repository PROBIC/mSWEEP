#ifndef MSWEEP_READ_PSEUDOALIGNMENT_HPP
#define MSWEEP_READ_PSEUDOALIGNMENT_HPP

#include <vector>
#include <memory>
#include <string>

#include "Sample.hpp"
#include "KallistoFiles.hpp"

void ReadPseudoalignment(KallistoFiles &kallisto_files, unsigned n_refs, std::vector<std::unique_ptr<Sample>> &batch, bool bootstrap_mode);
void ReadPseudoalignment(const std::string &tinfile1, const std::string &tinfile2, const std::string &themisto_mode, const bool bootstrap_mode, const unsigned n_refs, std::vector<std::unique_ptr<Sample>> &batch);
void ReadLikelihood(const bool bootstrap_mode, std::istream &likelihood_file, std::vector<std::unique_ptr<Sample>> &samples);

#endif

#ifndef MSWEEP_MSWEEP_HPP
#define MSWEEP_MSWEEP_HPP

#include <vector>
#include <memory>
#include <cstddef>
#include <string>

#include "DenseMatrix.hpp"

#include "Sample.hpp"
#include "Reference.hpp"
#include "Grouping.hpp"

void ReadGroupIndicators(const std::string &indicators_path, Reference *reference);

telescope::GroupedAlignment ReadPseudoalignments(const std::vector<std::string> &alignment_paths,
						 const std::string &themisto_merge_mode,
						 const Reference &reference);

seamat::DenseMatrix<double> ReadLikelihoodFromFile(const std::string &likelihood_path, const Reference &reference, std::ostream &log, std::vector<double> *log_ec_counts);

#endif

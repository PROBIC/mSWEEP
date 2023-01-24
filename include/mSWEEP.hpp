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
void WriteLikelihood(const seamat::DenseMatrix<double> &ll_mat, const std::vector<double> &log_ec_counts, const uint32_t n_groups, std::ostream &of);
void WriteLikelihoodBitSeq(const seamat::DenseMatrix<double> &ll_mat, const std::vector<double> &log_ec_counts, const uint32_t n_groups, std::ostream &of);
void WriteProbabilities(const seamat::DenseMatrix<double> &ec_probs, const std::vector<std::string> &cluster_indicators_to_string, std::ostream &of);
void WriteAbundances(const std::vector<double> &relative_abundances, const std::vector<std::string> &cluster_indicators_to_string, const size_t counts_total, std::ostream &of);
void WriteBootstrappedAbundances(const std::vector<double> &relative_abundances, const std::vector<std::vector<double>> &bootstrap_results, const std::vector<std::string> &cluster_indicators_to_string, const size_t counts_total, const uint16_t iters, std::ostream &of);

#endif

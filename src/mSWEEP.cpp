// mSWEEP: Estimate abundances of reference lineages in DNA sequencing reads.
//
// MIT License
//
// Copyright (c) 2023 Probabilistic Inference and Computational Biology group @ UH
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
#include "mSWEEP.hpp"

#include <sstream>
#include <exception>
#include <numeric>

#include "cxxio.hpp"

#include "version.h"

void ReadGroupIndicators(const std::string &indicators_path, Reference *reference) {
  cxxio::In indicators_file(indicators_path);
  reference->read_from_file(indicators_file.stream(), '\t'); // TODO: take delimiter as argument.
}

telescope::GroupedAlignment ReadPseudoalignments(const std::vector<std::string> &alignment_paths, const std::string &themisto_merge_mode, const Reference &reference) {
  size_t n_files = alignment_paths.size();
  std::vector<cxxio::In> infiles;
  infiles.reserve(n_files);
  std::vector<std::istream*> strands(n_files);
  if (n_files > 0) {
    for (size_t i = 0; i < n_files; ++i) {
      infiles.emplace_back(cxxio::In(alignment_paths[i]));
      strands[i] = &infiles[i].stream();
    }
  } else {
    strands.emplace_back(&std::cin);
  }
  return telescope::read::ThemistoGrouped(telescope::get_mode(themisto_merge_mode), reference.get_n_refs(), reference.get_group_indicators(0), strands);
}

seamat::DenseMatrix<double> ReadLikelihoodFromFile(const std::string &likelihood_path, const Reference &reference, std::ostream &log, std::vector<double> *log_ec_counts) {
  log << "  reading likelihoods from file" << '\n';
  if (reference.get_n_groupings() > 1) {
    throw std::runtime_error("Using more than one grouping with --read-likelihood is not yet implemented.");
  }
  cxxio::In infile(likelihood_path);

  // Have to read the likelihoods into a temporary because num_ecs is not known
  uint32_t n_groups = reference.get_grouping(0).get_n_groups();
  std::vector<std::vector<double>> likelihoods(n_groups, std::vector<double>());

  if (infile.stream().good()) {
    std::string newline;
    uint32_t line_nr = 0;
    while (std::getline(infile.stream(), newline)) {
      ++line_nr;
      std::string part;
      std::stringstream partition(newline);
      bool ec_count_col = true;
      uint32_t group_id = 0;
      while (std::getline(partition, part, '\t')) {
	if (ec_count_col) {
	  uint32_t ec_count = std::stol(part);
	  log_ec_counts->emplace_back(std::log(ec_count));
	  ec_count_col = false;
	} else {
	  likelihoods[group_id].emplace_back(std::stod(part));
	  ++group_id;
	}
      }
    }
  } else {
    throw std::runtime_error("Could not read from the likelihoods file.");
  }
  return seamat::DenseMatrix<double>(likelihoods);
}

void WriteLikelihood(const seamat::DenseMatrix<double> &ll_mat, const std::vector<double> &log_ec_counts, const uint32_t n_groups, std::ostream &of) {
  // Write likelihoods to a file
  if (of.good()) {
    for (uint32_t i = 0; i < log_ec_counts.size(); ++i){
      uint32_t ec_hit_count = std::round(std::exp(log_ec_counts[i]));
      of << ec_hit_count << '\t';
      for (uint32_t j = 0; j < n_groups; ++j) {
	of << ll_mat(j, i);
	of << (j == n_groups - 1 ? '\n' : '\t');
      }
    }
    of.flush();
  } else {
    throw std::runtime_error("Can't write to likelihoods file.");
  }
}

void WriteLikelihoodBitSeq(const seamat::DenseMatrix<double> &ll_mat, const std::vector<double> &log_ec_counts, const uint32_t n_groups, std::ostream &of) {
  // Write likelihoods to a file
  // *Note*: will write in BitSeq format!
  // Use Sample::write_likelihoods if tab-separated matrix format is needed.

  auto sum_of_exps = [](size_t accumulator, const double &val) {
    return accumulator + std::exp(val);
  };

  if (of.good()) {
    size_t counts_total = std::accumulate(log_ec_counts.begin(), log_ec_counts.end(), 0, sum_of_exps);
    of << "# Ntotal " << counts_total << '\n';
    of << "# Nmap " << counts_total << '\n';
    of << "# M " << n_groups << '\n';
    of << "# LOGFORMAT (probabilities saved on log scale.)" << '\n';
    of << "# r_name num_alignments (tr_id prob )^*{num_alignments}" << '\n';

    uint32_t read_id = 1;
    for (uint32_t i = 0; i < log_ec_counts.size(); ++i) {
      uint32_t ec_hit_count = std::round(std::exp(log_ec_counts[i]));
      for (uint32_t k = 0; k < ec_hit_count; ++k) {
	of << read_id << ' ';
	of << n_groups + 1 << ' ';
	for (uint32_t j = 0; j < n_groups; ++j) {
	  of << j + 1 << ' ' << ll_mat(j, i) << ' ';
	}
	of << 0 << ' ' << "-10000.00" << '\n';
	++read_id;
      }
    }
    of.flush();
  } else {
    throw std::runtime_error("Can't write to likelihoods file (bitseq format).");
  }
}

void WriteProbabilities(const seamat::DenseMatrix<double> &ec_probs, const std::vector<std::string> &cluster_indicators_to_string, std::ostream &of) {
  // Write the probability matrix to a file.
  if (of.good()) {
    of << "ec_id" << ',';
    size_t n_rows = ec_probs.get_rows();
    size_t n_cols = ec_probs.get_cols();
    for (uint32_t i = 0; i < n_rows; ++i) {
      of << cluster_indicators_to_string[i];
      of << (i < n_rows - 1 ? ',' : '\n');
    }
    for (uint32_t i = 0; i < n_cols; ++i) {
      of << i << ',';
      for (uint32_t j = 0; j < n_rows; ++j) {
	of << std::exp(ec_probs(j, i));
	of << (j < n_rows - 1 ? ',' : '\n');
      }
    }
    of << std::endl;
    of.flush();
  } else {
    throw std::runtime_error("Can't write to probs file.");
  }
}

void WriteAbundances(const std::vector<double> &relative_abundances, const std::vector<std::string> &cluster_indicators_to_string, const size_t counts_total, std::ostream &of) {
  // Write relative abundances to &of,
  if (of.good()) {
    of << "#mSWEEP_version:" << '\t' << MSWEEP_BUILD_VERSION << '\n';
    of << "#total_hits:" << '\t' << counts_total << '\n';
    of << "#c_id" << '\t' << "mean_theta" << '\n';
    for (size_t i = 0; i < relative_abundances.size(); ++i) {
      of << cluster_indicators_to_string[i] << '\t' << relative_abundances[i] << '\n';
    }
    of.flush();
  } else {
    throw std::runtime_error("Can't write to abundances file.");
  }
}

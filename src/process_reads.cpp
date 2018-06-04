#include "process_reads.hpp"

#include "likelihood.hpp"
#include "rcg.hpp"

void ProcessReads(const Reference &reference, const std::string &outfile, Sample &sample, OptimizerArgs args) {
  // Process pseudoalignments from kallisto.
  std::cerr << "Building log-likelihood array" << std::endl;

  const Matrix<double> &log_likelihood = likelihood_array_mat(sample, reference.grouping);

  std::cerr << "Estimating relative abundances" << std::endl;
  sample.ec_probs = rcg_optl_mat(log_likelihood, sample, args.alphas, args.tolerance, args.max_iters);

  sample.write_abundances(reference.group_names, outfile);  
  if (args.write_probs && !outfile.empty()) {
    sample.write_probabilities(reference.group_names, outfile);
  }
}

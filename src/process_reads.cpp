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

std::vector<double> ProcessReads2(const Reference &reference, Sample &sample, std::vector<long unsigned> ec_counts, OptimizerArgs args, unsigned iter) {
  // Process pseudoalignments but return the abundances rather than writing.
  std::cerr << "Building log-likelihood array" << std::endl;

  const Matrix<double> &log_likelihood = likelihood_array_mat(sample, reference.grouping);

  std::cerr << "Estimating relative abundances" << std::endl;
  sample.ec_probs = rcg_optl_mat(log_likelihood, sample.total_counts(), ec_counts, args.alphas, args.tolerance, args.max_iters);
  const std::vector<double> &abundances = sample.group_abundances();
  return abundances;
  //  sample.bootstrap_abundances.insert(std::make_pair(iter, abundances));
  //  std::cout << "inserting... " << abundances.size() << std::endl;
}

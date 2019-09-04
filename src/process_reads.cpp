#include <unordered_map>

#include "process_reads.hpp"
#include "likelihood.hpp"
#include "rcg.hpp"
#include "thread_pool.hpp"
#include "bootstrap.hpp"

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

void ProcessBatch(const Reference &reference, Arguments &args, std::vector<Sample> &bitfields) {
  // Don't launch extra threads if the batch is small
  args.nr_threads = (args.nr_threads > bitfields.size() ? bitfields.size() : args.nr_threads);
  ThreadPool pool(args.nr_threads);
  for (auto bitfield : bitfields) {
    std::string batch_outfile = (args.outfile.empty() ? args.outfile : args.outfile + "/" + bitfield.cell_name());
    pool.enqueue(&ProcessReads, reference, batch_outfile, bitfield, args.optimizer);
  }
}

void ProcessBootstrap(Reference &reference, Arguments &args, std::vector<Sample> &bitfields) {
  // Avoid launching extra threads if bootstrapping for only a few iterations
  args.nr_threads = (args.nr_threads > args.iters ? args.iters : args.nr_threads);
  ThreadPool pool(args.nr_threads);
  //  const std::unordered_map<std::string, std::vector<std::vector<double>>> &results = bootstrap_abundances(bitfields, reference, pool, args);
  const BootstrapResults &results = bootstrap_abundances(bitfields, reference, pool, args);

  for (auto kv : results.get()) {
    std::string outfile = (args.outfile.empty() || !args.batch_mode ? args.outfile : args.outfile + '/' + kv.first);
    pool.enqueue(&write_bootstrap, reference.group_names, kv.second.second, outfile, args.iters, kv.second.first);
  }
}

#include "process_reads.hpp"
#include "openmp_config.hpp"

#if defined(MSWEEP_OPENMP_SUPPORT) && (MSWEEP_OPENMP_SUPPORT) == 1
#include <omp.h>
#endif

#include <unordered_map>

#include "likelihood.hpp"
#include "rcg.hpp"
#include "thread_pool.hpp"
#include "bootstrap.hpp"
#include "zstr.hpp"

void ProcessReads(const Reference &reference, std::string outfile, Sample &sample, OptimizerArgs args) {
  // Process pseudoalignments from kallisto.
  std::cerr << "Building log-likelihood array" << std::endl;

#if defined(MSWEEP_OPENMP_SUPPORT) && (MSWEEP_OPENMP_SUPPORT) == 1
  omp_set_num_threads(args.nr_threads);
#endif
  const Matrix<double> &log_likelihood = likelihood_array_mat(sample, reference.grouping);

  std::cerr << "Estimating relative abundances" << std::endl;
  sample.ec_probs = rcg_optl_mat(log_likelihood, sample, args.alphas, args.tolerance, args.max_iters);

  sample.write_abundances(reference.group_names, outfile);  
  if (args.write_probs && !outfile.empty()) {
    std::unique_ptr<std::ostream> of;
    if (args.gzip_probs) {
      outfile += "_probs.csv.gz";
      of = std::unique_ptr<std::ostream>(new zstr::ofstream(outfile));
    } else {
      outfile += "_probs.csv";
      of = std::unique_ptr<std::ostream>(new std::ofstream(outfile));
    }
    sample.write_probabilities(reference.group_names, args.gzip_probs, (args.print_probs ? std::cout : *of));
  }
}

void ProcessBatch(const Reference &reference, Arguments &args, std::vector<Sample> &bitfields) {
  // Don't launch extra threads if the batch is small
  unsigned short nr_threads_estim = 1;
  if (args.optimizer.nr_threads > bitfields.size()) {
    nr_threads_estim = args.optimizer.nr_threads - bitfields.size(); // Let the estimation use more threads
    args.optimizer.nr_threads = bitfields.size();
  }
  ThreadPool pool(args.optimizer.nr_threads);
  args.optimizer.nr_threads = nr_threads_estim;
  for (auto bitfield : bitfields) {
    std::string batch_outfile = (args.outfile.empty() ? args.outfile : args.outfile + "/" + bitfield.cell_name());
    pool.enqueue(&ProcessReads, reference, batch_outfile, bitfield, args.optimizer);
  }
}

void ProcessBootstrap(Reference &reference, Arguments &args, std::vector<Sample> &bitfields) {
  // Avoid launching extra threads if bootstrapping for only a few iterations
  unsigned short nr_threads_estim = 1;
  if (args.optimizer.nr_threads > args.iters) {
    nr_threads_estim = args.optimizer.nr_threads - args.iters; // Let the estimation use more threads
    args.optimizer.nr_threads = args.iters;
  }
  ThreadPool pool(args.optimizer.nr_threads);
  args.optimizer.nr_threads = nr_threads_estim;
  //  const std::unordered_map<std::string, std::vector<std::vector<double>>> &results = bootstrap_abundances(bitfields, reference, pool, args);
  const BootstrapResults &results = bootstrap_abundances(bitfields, reference, pool, args);

  for (auto kv : results.get()) {
    std::string outfile = (args.outfile.empty() || !args.batch_mode ? args.outfile : args.outfile + '/' + kv.first);
    pool.enqueue(&write_bootstrap, reference.group_names, kv.second.second, outfile, args.iters, kv.second.first);
  }
}

#include "process_reads.hpp"
#include "openmp_config.hpp"

#include <unordered_map>

#include "likelihood.hpp"
#include "rcg.hpp"
#include "thread_pool.hpp"
#include "zstr.hpp"

void ProcessReads(const Reference &reference, std::string outfile, Sample &sample, OptimizerArgs args) {
  // Process pseudoalignments from kallisto.
  std::cerr << "Building log-likelihood array" << std::endl;

  likelihood_array_mat(sample, reference.grouping);

  std::cerr << "Estimating relative abundances" << std::endl;
  sample.ec_probs = rcg_optl_mat(sample.ll_mat, sample, args.alphas, args.tolerance, args.max_iters);

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
  for (auto bitfield : bitfields) {
    std::string batch_outfile = (args.outfile.empty() ? args.outfile : args.outfile + "/" + bitfield.cell_name());
    ProcessReads(reference, batch_outfile, bitfield, args.optimizer);
  }
}

void ProcessBootstrap(Reference &reference, Arguments &args, std::vector<Sample> &bitfields) {
  for (auto bitfield : bitfields) {
    BootstrapSample* bs = static_cast<BootstrapSample*>(&bitfield);
    bs->BootstrapAbundances(reference, args);
    bs->WriteBootstrap(reference.group_names, args.outfile, args.iters, args.batch_mode);    
  }
}

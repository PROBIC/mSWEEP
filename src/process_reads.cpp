#include "process_reads.hpp"

#include "rcg.hpp"
#include "bxzstr.hpp"

void ProcessReads(const Grouping &grouping, std::string outfile, Sample &sample, OptimizerArgs args) {
  // Process pseudoalignments.
  sample.ec_probs = rcg_optl_mat(sample.ll_mat, sample, args.alphas, args.tolerance, args.max_iters);

  sample.write_abundances(grouping.names, outfile);
  if (args.write_probs && !outfile.empty()) {
    std::unique_ptr<std::ostream> of;
    if (args.gzip_probs) {
      outfile += "_probs.csv.gz";
      of = std::unique_ptr<std::ostream>(new bxz::ofstream(outfile));
    } else {
      outfile += "_probs.csv";
      of = std::unique_ptr<std::ostream>(new std::ofstream(outfile));
    }
    sample.write_probabilities(grouping.names, args.gzip_probs, (args.print_probs ? std::cout : *of));
  }
}

void ProcessBatch(const Grouping &grouping, Arguments &args, std::vector<std::unique_ptr<Sample>> &bitfields) {
  for (uint32_t i = 0; i < bitfields.size(); ++i) {
    std::string batch_outfile = (args.outfile.empty() ? args.outfile : args.outfile + "/" + bitfields[i]->cell_name());
    ProcessReads(grouping, batch_outfile, *bitfields[i], args.optimizer);
  }
}

void ProcessBootstrap(const Grouping &grouping, Arguments &args, std::vector<std::unique_ptr<Sample>> &bitfields) {
  for (uint32_t i = 0; i < bitfields.size(); ++i) {
    BootstrapSample* bs = static_cast<BootstrapSample*>(&(*bitfields[i]));
    bs->BootstrapAbundances(grouping, args);
    bs->WriteBootstrap(grouping.names, args.outfile, args.iters, args.batch_mode);    
  }
}

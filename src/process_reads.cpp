#include "process_reads.hpp"

#include <string>

#include "rcg.hpp"
#include "bxzstr.hpp"

void ProcessReads(const Grouping &grouping, const Arguments &args, std::vector<std::unique_ptr<Sample>> &samples) {
  for (uint32_t i = 0; i < samples.size(); ++i) {
    // Process pseudoalignments.
    samples[i]->ec_probs = rcg_optl_mat(samples[i]->ll_mat, (*samples[i]), args.optimizer.alphas, args.optimizer.tolerance, args.optimizer.max_iters);

    std::string outfile(args.outfile);
    if (samples.size() > 1) {
      // Legacy kallisto batch mode support in outfile names.
      outfile = (args.outfile.empty() ? args.outfile : args.outfile + "/" + samples[i]->cell_name());
    }

    samples[i]->write_abundances(grouping.get_names(), outfile);
    if (args.optimizer.write_probs && !outfile.empty()) {
      std::unique_ptr<std::ostream> of;
      if (args.optimizer.gzip_probs) {
	outfile += "_probs.csv.gz";
	of = std::unique_ptr<std::ostream>(new bxz::ofstream(outfile));
      } else {
	outfile += "_probs.csv";
	of = std::unique_ptr<std::ostream>(new std::ofstream(outfile));
      }
      samples[i]->write_probabilities(grouping.get_names(), (args.optimizer.print_probs ? std::cout : *of));
    }
  }
}

void ProcessBootstrap(const Grouping &grouping, Arguments &args, std::vector<std::unique_ptr<Sample>> &bitfields) {
  for (uint32_t i = 0; i < bitfields.size(); ++i) {
    BootstrapSample* bs = static_cast<BootstrapSample*>(&(*bitfields[i]));
    bs->BootstrapAbundances(grouping, args);
    bs->WriteBootstrap(grouping.get_names(), args.outfile, args.iters, args.batch_mode);    
  }
}

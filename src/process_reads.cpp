#include "process_reads.hpp"

#include "rcg.hpp"
#include "bxzstr.hpp"

void ProcessReads(const Reference &reference, std::string outfile, Sample &sample, OptimizerArgs args) {
  // Process pseudoalignments from kallisto.
  std::cerr << "Building log-likelihood array" << std::endl;

  sample.CalcLikelihood(reference.grouping);

  if (args.write_likelihood || args.write_likelihood_bitseq) {
    std::cerr << "Writing likelihood matrix" << std::endl;
    if (args.write_likelihood) {
      sample.write_likelihood(args.gzip_probs, reference.grouping.n_groups, outfile);
    }
    if (args.write_likelihood_bitseq) {
      sample.write_likelihood_bitseq(args.gzip_probs, reference.grouping.n_groups, outfile);
    }
  }

  if (args.no_fit_model) {
    std::cerr << "Skipping relative abundance estimation (--no-fit-model toggled)" << std::endl;
  } else {
    std::cerr << "Estimating relative abundances" << std::endl;
    sample.ec_probs = rcg_optl_mat(sample.ll_mat, sample, args.alphas, args.tolerance, args.max_iters);

    sample.write_abundances(reference.group_names, outfile);  
    if (args.write_probs && !outfile.empty()) {
      std::unique_ptr<std::ostream> of;
      if (args.gzip_probs) {
	outfile += "_probs.csv.gz";
	of = std::unique_ptr<std::ostream>(new bxz::ofstream(outfile));
      } else {
	outfile += "_probs.csv";
	of = std::unique_ptr<std::ostream>(new std::ofstream(outfile));
      }
      sample.write_probabilities(reference.group_names, args.gzip_probs, (args.print_probs ? std::cout : *of));
    }
  }
}

void ProcessBatch(const Reference &reference, Arguments &args, std::vector<std::unique_ptr<Sample>> &bitfields) {
  for (uint32_t i = 0; i < bitfields.size(); ++i) {
    std::string batch_outfile = (args.outfile.empty() ? args.outfile : args.outfile + "/" + bitfields[i]->cell_name());
    ProcessReads(reference, batch_outfile, *bitfields[i], args.optimizer);
  }
}

void ProcessBootstrap(Reference &reference, Arguments &args, std::vector<std::unique_ptr<Sample>> &bitfields) {
  for (uint32_t i = 0; i < bitfields.size(); ++i) {
    BootstrapSample* bs = static_cast<BootstrapSample*>(&(*bitfields[i]));
    bs->BootstrapAbundances(reference, args);
    bs->WriteBootstrap(reference.group_names, args.outfile, args.iters, args.batch_mode);    
  }
}

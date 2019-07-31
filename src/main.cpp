#include <iostream>
#include <vector>
#include <exception>

#include "parse_arguments.hpp"
#include "read_bitfield.hpp"
#include "process_reads.hpp"
#include "thread_pool.hpp"
#include "Sample.hpp"
#include "Reference.hpp"
#include "version.h"

int main (int argc, char *argv[]) {
  std::cerr << "mSWEEP-" << _BUILD_VERSION << " abundance estimation" << std::endl;
  Arguments args;
  try {
    ParseArguments(argc, argv, args);
  }
  catch (std::runtime_error &e) {
    std::cerr << "Error in parsing arguments:\n  "
	      << e.what()
	      << "\nexiting" << std::endl;
    return 1;
  }
  catch (std::invalid_argument &e) {
    std::cerr << e.what() << std::endl;
    PrintHelpMessage();
    return 0;
  }

  std::vector<Sample> bitfields;
  Reference reference;
  try {
    std::cerr << "Reading the input files" << '\n';
    std::cerr << "  reading group indicators" << '\n';
    ReadClusterIndicators(args.indicators_file, reference);
    std::cerr << "  read " << reference.n_refs << " group indicators" << std::endl;

    // Check that the number of reference sequences matches in the grouping and the alignment.
    VerifyGrouping(args.kallisto_files[0], reference.n_refs);

    std::cerr << "  reading pseudoalignments" << '\n';
    ReadBitfield(args.kallisto_files, reference.n_refs, bitfields);
    std::cerr << "  read " << (args.batch_mode ? bitfields.size() : bitfields[0].num_ecs()) << (args.batch_mode ? " samples from the batch" : " unique alignments") << std::endl;
  } catch (std::runtime_error &e) {
    std::cerr << "Reading pseudoalignments failed:\n  ";
    std::cerr << e.what();
    std::cerr << "\nexiting" << std::endl;
    return 1;
  }

  // Calculate the beta-binomial parameters for the grouping
  reference.calculate_bb_parameters(args.params);

  // Initialize the prior counts on the groups
  args.optimizer.alphas = std::vector<double>(reference.grouping.n_groups, 1.0);

  if (!args.batch_mode && args.iters == 1) {
    ProcessReads(reference, args.outfile, bitfields[0], args.optimizer);
  } else if (args.iters == 1) {
    ProcessBatch(reference, args, bitfields);
  } else {
    ProcessBootstrap(reference, args, bitfields);
  }

  return 0;
}

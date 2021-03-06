#include <iostream>
#include <vector>
#include <exception>
#include <memory>

#include "bxzstr.hpp"
#include "file.hpp"

#include "parse_arguments.hpp"
#include "read_bitfield.hpp"
#include "process_reads.hpp"
#include "Sample.hpp"
#include "Reference.hpp"
#include "version.h"
#include "openmp_config.hpp"

#if defined(MSWEEP_OPENMP_SUPPORT) && (MSWEEP_OPENMP_SUPPORT) == 1
#include <omp.h>
#endif

int main (int argc, char *argv[]) {
  std::cerr << "mSWEEP-" << MSWEEP_BUILD_VERSION << " abundance estimation" << std::endl;
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

#if defined(MSWEEP_OPENMP_SUPPORT) && (MSWEEP_OPENMP_SUPPORT) == 1
  omp_set_num_threads(args.optimizer.nr_threads);
#endif

  std::vector<std::unique_ptr<Sample>> bitfields;
  Reference reference;
  try {
    std::cerr << "Reading the input files" << '\n';
    std::cerr << "  reading group indicators" << '\n';
    if (args.fasta_file.empty()) {
      File::In indicators_file(args.indicators_file);
      ReadClusterIndicators(indicators_file.stream(), reference);
    } else {
      File::In groups_file(args.groups_list_file);
      File::In fasta_file(args.fasta_file);
      MatchClusterIndicators(args.groups_list_delimiter, groups_file.stream(), fasta_file.stream(), reference);
    }
    if (reference.n_refs == 0) {
      throw std::runtime_error("The grouping contains 0 reference sequences");
    }
    std::cerr << "  read " << reference.n_refs << " group indicators" << std::endl;

    std::cerr << "  reading pseudoalignments" << '\n';
    if (!args.themisto_mode) {
      // Check that the number of reference sequences matches in the grouping and the alignment.
      VerifyGrouping(reference.n_refs, *args.infiles.run_info);
      ReadBitfield(args.infiles, reference.n_refs, bitfields, reference, args.bootstrap_mode);
    } else {
      if (!args.themisto_index_path.empty()) {
	File::In themisto_index(args.themisto_index_path + "/coloring-names.txt");
	VerifyThemistoGrouping(reference.n_refs, themisto_index.stream());
      }
      ReadBitfield(args.tinfile1, args.tinfile2, args.themisto_merge_mode, args.bootstrap_mode, reference.n_refs, bitfields);
    }

    std::cerr << "  read " << (args.batch_mode ? bitfields.size() : bitfields[0]->num_ecs()) << (args.batch_mode ? " samples from the batch" : " unique alignments") << std::endl;
  } catch (std::runtime_error &e) {
    std::cerr << "Reading the input files failed:\n  ";
    std::cerr << e.what();
    std::cerr << "\nexiting" << std::endl;
    return 1;
  }

  // Calculate the beta-binomial parameters for the grouping
  reference.calculate_bb_parameters(args.params);

  // Initialize the prior counts on the groups
  args.optimizer.alphas = std::vector<double>(reference.grouping.n_groups, 1.0);

  // Process the reads accordingly
  switch(args.run_mode()) {
  case 0: ProcessReads(reference, args.outfile, *bitfields[0], args.optimizer); break;
  case 1: ProcessBatch(reference, args, bitfields); break;
  case 2: ProcessBootstrap(reference, args, bitfields); break;
  case 3: ProcessBootstrap(reference, args, bitfields); break; // Same function for batch and single files
  }

  return 0;
}

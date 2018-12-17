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

  std::cerr << "Reading the input files" << '\n';
  bool batch_mode = args.infile.empty();
  std::vector<std::string> kallisto_files((batch_mode ? 4 : 3));
  if (batch_mode) {
    kallisto_files[0] = args.batch_infile + "/run_info.json";
    kallisto_files[1] = args.batch_infile + "/matrix.ec";
    kallisto_files[2] = args.batch_infile + "/matrix.tsv";
    kallisto_files[3] = args.batch_infile + "/matrix.cells";
  } else {
    kallisto_files[0] = args.infile + "/run_info.json";
    kallisto_files[1] = args.infile + "/pseudoalignments.ec";
    kallisto_files[2] = args.infile + "/pseudoalignments.tsv";
  }

  std::vector<Sample> bitfields;
  Reference reference;
  try {
    std::cerr << "  reading group indicators" << '\n';
    ReadClusterIndicators(args.indicators_file, reference);
    std::cerr << "  read " << reference.n_refs << " group indicators" << std::endl;

    // Check that the number of reference sequences matches in the grouping and the alignment.
    VerifyGrouping(kallisto_files[0], reference.n_refs);

    std::cerr << "  reading pseudoalignments" << '\n';
    ReadBitfield(kallisto_files, reference.n_refs, bitfields);
    std::cerr << "  read " << (batch_mode ? bitfields.size() : bitfields[0].num_ecs()) << (batch_mode ? " samples from the batch" : " unique alignments") << std::endl;
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

  if (!batch_mode) {
    ProcessReads(reference, args.outfile, bitfields[0], args.optimizer);
  } else {
    // Don't launch extra threads if the batch is small
    args.nr_threads = (args.nr_threads > bitfields.size() ? bitfields.size() : args.nr_threads);
    ThreadPool pool(args.nr_threads);
    for (auto bitfield : bitfields) {
      std::string batch_outfile = (args.outfile.empty() ? args.outfile : args.outfile + "/" + bitfield.cell_name());
      pool.enqueue(&ProcessReads, reference, batch_outfile, bitfield, args.optimizer);
    }
  }
  
  return 0;
}

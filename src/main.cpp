#include <iostream>
#include <vector>
#include <exception>
#include <random>
#include <future>

#include "parse_arguments.hpp"
#include "read_bitfield.hpp"
#include "process_reads.hpp"
#include "thread_pool.hpp"
#include "Sample.hpp"
#include "Reference.hpp"

void write_bootstrap(const std::vector<std::string> &cluster_indicators_to_string, const std::vector<std::vector<double>> &abundances, std::string &outfile, unsigned iters) {
  // Write relative abundances to a file,
  // outputs to std::cout if outfile is empty.
  std::streambuf *buf;
  std::ofstream of;
  if (outfile.empty()) {
    buf = std::cout.rdbuf();
  } else {
    outfile += "_abundances.txt";
    of.open(outfile);
    buf = of.rdbuf();
  }
  std::ostream out(buf);
  out << "#c_id" << '\t' << "abundances" << '\t' << "bootstrap_abundances" << '\n';

  for (size_t i = 0; i < cluster_indicators_to_string.size(); ++i) {
    out << cluster_indicators_to_string[i] << '\t';
    for (unsigned j = 0; j < iters; ++j) {
      out << abundances.at(j).at(i) << (j == iters - 1 ? '\n' : '\t');
    }
  }
  out << std::endl;
  if (!outfile.empty()) {
    of.close();
  }
}

int main (int argc, char *argv[]) {
  std::cerr << "mSWEEP relative abundance estimation" << std::endl;
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

  if (!batch_mode && args.iters == 1) {
    ProcessReads(reference, args.outfile, bitfields[0], args.optimizer);
  } else if (args.iters == 1) {
    // Don't launch extra threads if the batch is small
    args.nr_threads = (args.nr_threads > bitfields.size() ? bitfields.size() : args.nr_threads);
    ThreadPool pool(args.nr_threads);
    for (auto bitfield : bitfields) {
      std::string batch_outfile = (args.outfile.empty() ? args.outfile : args.outfile + "/" + bitfield.cell_name());
      pool.enqueue(&ProcessReads, reference, batch_outfile, bitfield, args.optimizer);
    }
  } else {
    std::unordered_map<std::string, std::vector<std::vector<double>>> results;
    // There's probably a more elegant way to implement this
    args.nr_threads = (args.nr_threads > args.iters ? args.iters : args.nr_threads);
    ThreadPool pool(args.nr_threads);
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::cerr << "Running estimation with " << args.iters << " bootstrap iterations" << '\n';
    for (auto bitfield : bitfields) {
      // Store results in this
      std::vector<std::future<std::vector<double>>> abus;
      // Init the bootstrap variables
      bitfield.init_bootstrap();
      for (unsigned i = 0; i < args.iters; ++i) {
	// Run the estimation multiple times without writing anything
	abus.emplace_back(pool.enqueue(&ProcessReads2, reference, bitfield, bitfield.ec_counts, args.optimizer, i));
	// Resample the pseudoalignment counts (here because we want to include the original)
	bitfield.resample_counts(gen);
      }
      std::string name = (batch_mode ? bitfield.cell_name() : "0");
      results.insert(std::make_pair(name, std::vector<std::vector<double>>()));
      for (unsigned i = 0; i < args.iters; ++i) {
	results.at(name).emplace_back(abus.at(i).get());
      }
    }
    for (auto kv : results) {
      std::string outfile = (args.outfile.empty() || !batch_mode ? args.outfile : args.outfile + '/' + kv.first);
      pool.enqueue(&write_bootstrap, reference.group_names, kv.second, outfile, args.iters);
    }
  }

  return 0;
}

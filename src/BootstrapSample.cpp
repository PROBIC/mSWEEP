#include "Sample.hpp"

#include "likelihood.hpp"
#include "rcg.hpp"
#include "version.h"

#include "bxzstr.hpp"

void BootstrapSample::InitBootstrap(const Grouping &grouping, const double bb_constants[2], const std::vector<uint32_t> &group_indicators) {
  ec_distribution = std::discrete_distribution<uint32_t>(pseudos.ec_counts.begin(), pseudos.ec_counts.end());
  CalcLikelihood(grouping, bb_constants, group_indicators);

  // Clear the abundances in case we're estimating the same sample again.
  this->relative_abundances = std::vector<std::vector<double>>();
}

void BootstrapSample::ResampleCounts(const uint32_t how_many, std::mt19937_64 &generator) {
  std::vector<uint32_t> tmp_counts(num_ecs());
  for (uint32_t i = 0; i < how_many; ++i) {
    uint32_t ec_id = ec_distribution(generator);
    tmp_counts[ec_id] += 1;
  }
#pragma omp parallel for schedule(static)
  for (uint32_t i = 0; i < num_ecs(); ++i) {
    log_ec_counts[i] = std::log(tmp_counts[i]);
  }
  counts_total = how_many;
}

void BootstrapSample::BootstrapIter(const std::vector<double> &alpha0, const double tolerance, const uint16_t max_iters) {
  // Process pseudoalignments but return the abundances rather than writing.
  ec_probs = rcg_optl_mat(ll_mat, *this, alpha0, tolerance, max_iters);
  this->relative_abundances.emplace_back(group_abundances());
}

void BootstrapSample::BootstrapAbundances(const Grouping &grouping, const std::vector<uint32_t> &group_indicators, const Arguments &args) {
  std::mt19937_64 gen;
  if (args.seed == -1) {
    std::random_device rd;
    gen = std::mt19937_64(rd());
  } else {
    gen = std::mt19937_64(args.seed);
  }
  std::cerr << "Running estimation with " << args.iters << " bootstrap iterations" << '\n';
  // Which sample are we processing?
  std::string name = (args.batch_mode ? cell_name() : "0");
  std::cout << "Processing " << (args.batch_mode ? name : "the sample") << std::endl;
  // Init the bootstrap variables
  InitBootstrap(grouping, args.optimizer.bb_constants, group_indicators);

  // Store the original values
  std::vector<double> og_log_ec_counts(this->log_ec_counts);
  uint32_t og_counts_total = this->counts_total;

  for (uint16_t i = 0; i <= args.iters; ++i) {
    if (i > 0) {
      std::cout << "Bootstrap" << " iter " << i << "/" << args.iters << std::endl;
    } else {
      std::cerr << "Estimating relative abundances without bootstrapping" << std::endl;
    }
    BootstrapIter(args.optimizer.alphas, args.optimizer.tolerance, args.optimizer.max_iters);

    if (i == 0) {
      if (args.optimizer.write_probs && !args.outfile.empty()) {
	std::string outfile = args.outfile;
	std::unique_ptr<std::ostream> of;
	if (args.optimizer.gzip_probs) {
	  outfile += "_probs.csv.gz";
	  of = std::unique_ptr<std::ostream>(new bxz::ofstream(outfile));
	} else {
	  outfile += "_probs.csv";
	  of = std::unique_ptr<std::ostream>(new std::ofstream(outfile));
	}
	write_probabilities(grouping.names, args.optimizer.gzip_probs, (args.optimizer.print_probs ? std::cout : *of));
      }
    }
    // Resample the pseudoalignment counts (here because we want to include the original)
    ResampleCounts((args.bootstrap_count == 0 ? counts_total : args.bootstrap_count), gen);
  }

  // Restore original values
  this->log_ec_counts = og_log_ec_counts;
  this->counts_total = og_counts_total;
}


void BootstrapSample::WriteBootstrap(const std::vector<std::string> &cluster_indicators_to_string, std::string &outfile, const unsigned iters, const bool batch_mode) const {
  // Write relative abundances to a file,
  // outputs to std::cout if outfile is empty.
  outfile = (outfile.empty() || !batch_mode ? outfile : outfile + '/' + cell_name());
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
  out << "#mSWEEP_version:" << '\t' << MSWEEP_BUILD_VERSION << '\n';
  out << "#total_hits:" << '\t' << counts_total << '\n';
  out << "#bootstrap_iters:" << '\t' << iters << '\n';
  out << "#c_id" << '\t' << "mean_theta" << '\t' << "bootstrap_mean_thetas" << '\n';

  for (size_t i = 0; i < cluster_indicators_to_string.size(); ++i) {
    out << cluster_indicators_to_string[i] << '\t';
    for (unsigned j = 0; j <= iters; ++j) {
      out << relative_abundances[j][i] << (j == iters ? '\n' : '\t');
    }
  }
  if (!outfile.empty()) {
    of.close();
  }
}

void BootstrapSample::read_themisto(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &strands) {
  ReadThemisto(mode, n_refs, strands, &pseudos);
  process_aln();
}

void BootstrapSample::read_kallisto(const uint32_t n_refs, std::istream &ec_file, std::istream &tsv_file) {
  ReadKallisto(n_refs, ec_file, tsv_file, &pseudos);
  process_aln();
}


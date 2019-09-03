#ifndef BOOTSTRAP_CONFINTS_H
#define BOOTSTRAP_CONFINTS_H

#include <string>
#include <vector>
#include <future>
#include <random>

#include "Sample.hpp"
#include "Reference.hpp"
#include "thread_pool.hpp"
#include "parse_arguments.hpp"

void write_bootstrap(const std::vector<std::string> &cluster_indicators_to_string, const std::vector<std::vector<double>> &abundances, std::string &outfile, unsigned iters, unsigned counts_total);
std::vector<std::vector<double>> bootstrap_abundances(const Arguments &args, std::mt19937_64 &gen, Reference &reference, Sample &sample, ThreadPool &pool);

#endif

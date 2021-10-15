#ifndef MSWEEP_PROCESS_READS_HPP
#define MSWEEP_PROCESS_READS_HPP

#include <vector>
#include <memory>

#include "parse_arguments.hpp"
#include "Reference.hpp"
#include "Sample.hpp"

void ProcessReads(const Grouping &grouping, const Arguments &args, std::vector<std::unique_ptr<Sample>> &samples);
void ProcessBootstrap(const Grouping &grouping, Arguments &args, std::vector<std::unique_ptr<Sample>> &bitfields);

#endif

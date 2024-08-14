// mSWEEP: Estimate abundances of reference lineages in DNA sequencing reads.
//
// MIT License
//
// Copyright (c) 2023 Probabilistic Inference and Computational Biology group @ UH
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
#ifndef MSWEEP_ALIGNMENT_HPP
#define MSWEEP_ALIGNMENT_HPP

#include "bm64.h"
#include "unpack.hpp"

#include <fstream>
#include <sstream>
#include <iostream>

namespace mSWEEP {
class Alignment {
private:
    bm::bvector<> bits;
    size_t n_targets;
    size_t n_queries;
    std::vector<size_t> group_indicators;
    size_t n_groups;

public:
    Alignment(size_t _n_targets) {
	this->n_targets = _n_targets;
    }

    void ReadPlaintextLine(const size_t n_targets, std::string &line, bm::bvector<>::bulk_insert_iterator &it) {
	// telescope::ReadPlaintextLine
	//
	// Reads a line in a plaintext alignment file from Themisto
	// (https://github.com/algbio/themisto) into the bm::bvector<> `it`
	// inserts to.
	//
	// Input:
	//   `n_targets`: number of pseudoalignment targets (reference
	//                sequences). It's not possible to infer this from the Themisto
	//                file format so has to be provided separately.
	//   `line`: the line from the alignment file to read in.
	//   `it`: insert iterator to the bm::bvector<> variable for storing the alignment.
	//
	std::string part;
	std::stringstream partition(line);

	// First column is read id (0-based indexing).
	std::getline(partition, part, ' ');
	size_t read_id = std::stoul(part);

	// Next columns contain the target sequence id (0-based indexing).
	while (std::getline(partition, part, ' ')) {
	    *it = read_id*n_targets + std::stoul(part); // set bit `n_reads*n_refs + std::stoul(part)` as true
	}
    }

    size_t ReadPlaintextAlignment(const size_t n_targets, std::string &line, std::istream *stream, bm::bvector<> *ec_configs) {
	// telescope::ReadPlaintextAlignment
	//
	// Reads a plaintext alignment file from Themisto
	// (https://github.com/algbio/themisto) into `*ec_configs` and
	// return the number of reads in the file (both unaligned and
	// aligned).
	//
	// Input:
	//   `n_targets`: number of pseudoalignment targets (reference
	//                sequences). It's not possible to infer this from the Themisto
	//                file format so has to be provided separately.
	//   `line`: read each line into this variable.
	//     NOTE: the contents of the *first* line should already be stored in this variable.
	//   `stream`: pointer to an istream opened on the pseudoalignment file.
	//   `ec_configs`: pointer to the output variable that will contain the alignment.
	// Output:
	//   `n_reads`: total number of reads in the pseudoalignment (unaligned + aligned).
	//
	bm::bvector<>::bulk_insert_iterator it(*ec_configs); // Bulk insert iterator buffers the insertions

	size_t n_reads = 1;
	try {
	    // Contents of the first line is already stored in `line`
	    ReadPlaintextLine(n_targets, line, it);

	    size_t compress_interval = 1000000;
	    while (std::getline(*stream, line)) {
		// Insert each line into the alignment
		ReadPlaintextLine(n_targets, line, it);
		++n_reads;
		if (n_reads % compress_interval == 0) {
		    ec_configs->optimize();
		}
	    }
	} catch (const std::exception &e) {
	    std::string msg(e.what());
	    if (msg.find("stoul") != std::string::npos) {
		throw std::runtime_error("File format not supported on line " + std::to_string(n_reads) + " with content: " + line);
	    } else {
		throw std::runtime_error("Could not parse line " + std::to_string(n_reads) + " with content: " + line);
	    }
	}
	return n_reads;
    }


    void read(const std::string &merge_mode, std::vector<std::istream*> &strands) {
	for (size_t i = 0; i < strands.size(); ++i) {
	    std::string line;
	    std::getline(*strands[i], line); // Read the first line to check the format
	    size_t n_reads;
	    bm::bvector<> strand_alignment;
	    if (line.find(',') != std::string::npos) {
		// First line contains a ','; stream could be in the compact format.
		size_t n_refs;
		alignment_writer::ReadHeader(line, &n_reads, &n_refs);
		if (n_refs > this->n_targets) {
		    throw std::runtime_error("Pseudoalignment file has more target sequences than expected.");
		} else if (this->n_targets < n_refs) {
		    throw std::runtime_error("Pseudoalignment file has less target sequences than expected.");
		}
		// Size is given on the header line.
		strand_alignment.resize(n_reads*n_refs);
		alignment_writer::UnpackData(strands[i], strand_alignment);
	    } else {
		// Stream could be in the plaintext format.
		// Size is unknown.
		strand_alignment.set_new_blocks_strat(bm::BM_GAP);
		n_reads = ReadPlaintextAlignment(n_targets, line, strands[i], &strand_alignment);
	    }
	    this->n_queries = n_reads;

	    if (i == 0) {
		this->bits = std::move(strand_alignment);
	    } else {
		if (merge_mode == "intersection") {
		    this->bits.bit_and(strand_alignment);
		} else if (merge_mode == "union") {
		    this->bits.bit_or(strand_alignment);
		} else {
		    throw std::runtime_error("Unrecognized option `" + merge_mode + "` for --themisto-mode");
		}
	    }
	}
    }

    void collapse() {
	
    }

    size_t n_ecs() const { return this->n_queries; };
    size_t n_reads() const { return this->n_queries; };
    size_t reads_in_ec(const size_t i=0) const { return 1; };

    std::vector<size_t> operator()(const size_t row) const {
	std::vector<size_t> ret(this->n_groups, 0);
	for (size_t i = 0; i < this->n_targets; ++i) {
	    ret[this->group_indicators[i]] += this->bits[row*n_targets + i];
	}
	return ret;
    }

    std::vector<std::vector<uint32_t>> get_aligned_reads() const {
	std::vector<std::vector<uint32_t>> ret(this->n_queries, std::vector<uint32_t>(1, 0));
	for (size_t i = 0; i < n_queries; ++i) {
	    ret[i][0] = i;
	}
	return ret;
    }

    template <typename T>
    void add_groups(const std::vector<T> &grouping) {
	size_t _n_groups = 0;
	this->group_indicators = std::vector<size_t>(grouping.size());
	for (size_t i = 0; i < grouping.size(); ++i) {
	    group_indicators[i] = grouping[i];
	    _n_groups = (n_groups > grouping[i] ? n_groups : grouping[i]);
	}
	this->n_groups = _n_groups + 1;
    }

};
}

#endif

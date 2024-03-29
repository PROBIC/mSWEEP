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
#include "OutfileDesignator.hpp"

#include <exception>

namespace mSWEEP {
void OutfileDesignator::open(std::string &filename) {
  if (this->compress) {
    filename += this->extension;
    of.open_compressed(filename, this->type, this->compression_level);
  } else {
    of.open(filename);
  }
}

OutfileDesignator::OutfileDesignator(std::string _prefix, size_t _n_groupings, std::string _compress, int _level) {
  this->printing = _prefix.empty();
  this->prefix = _prefix;
  this->n_groupings = _n_groupings;

  this->compress = (_compress != "plaintext");
  if (this->compress) {
    this->compression_level = _level;
    if (_compress == "z") {
      this->type = bxz::z;
      this->extension = ".gz";
    } else if (_compress == "bz2") {
      this->type = bxz::bz2;
      this->extension = ".bz2";
    } else if (_compress == "lzma") {
      this->type = bxz::lzma;
      this->extension = ".xz";
    } else if (_compress == "zstd") {
      this->type = bxz::zstd;
      this->extension = ".zst";
    } else {
      throw std::invalid_argument("unsupported compression type " + _compress);
    }
  }

  if (this->n_groupings > 1) {
    this->current_grouping = 0;
    this->prefix += "_";
    this->prefix += std::to_string(this->current_grouping);
  }
}

std::ostream* OutfileDesignator::likelihoods(const std::string &format) {
  std::string ll_outfile = this->prefix;
  ll_outfile += (format == "bitseq" ? "_bitseq" : "");
  ll_outfile += "_likelihoods.tsv";

  this->open(ll_outfile);
  return &this->of.stream();
}

std::ostream* OutfileDesignator::bin(const std::string &name) {
  std::string bin_outfile;
  if (this->prefix.find('/') != std::string::npos) {
    // If the outfile location is in another folder then get the path
    bin_outfile = this->prefix;
    bin_outfile.erase(bin_outfile.rfind("/"), bin_outfile.size());
  } else {
    // If not in a folder write into the current directory.
    bin_outfile = ".";
  }

  bin_outfile += '/' + name + ".bin";
  this->open(bin_outfile);
  return &this->of.stream();
}

std::ostream* OutfileDesignator::probs() {
  std::string probs_outfile = this->prefix;
  probs_outfile += "_probs.tsv";

  this->open(probs_outfile);
  return &this->of.stream();
}

std::ostream* OutfileDesignator::abundances() {
  if (!printing) {
    std::string abundances_outfile = this->prefix + "_abundances.txt";
    this->of.open(abundances_outfile); // Ignore request to compress
  } else {
    if (&this->of.stream() != &std::cout) {
      this->of.close();
    }
  }
  return &this->of.stream();
}

void OutfileDesignator::next_grouping() {
  ++this->current_grouping;
  if (!printing) {
    this->prefix.erase(this->prefix.rfind("_"), this->prefix.size());
    this->prefix += "_";
    this->prefix += std::to_string(this->current_grouping);
  }
}
}

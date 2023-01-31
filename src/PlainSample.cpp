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
#include "Sample.hpp"

#include <exception>

#include "mSWEEP_version.h"

namespace mSWEEP {
void PlainSample::write_abundances(const std::vector<std::string> &group_names, std::ostream *of) const {
  // Write relative abundances to &of,
  if (of->good()) {
    (*of) << "#mSWEEP_version:" << '\t' << MSWEEP_BUILD_VERSION << '\n';
    (*of) << "#total_hits:" << '\t' << this->get_counts_total() << '\n';
    (*of) << "#c_id" << '\t' << "mean_theta" << '\n';
    for (size_t i = 0; i < this->relative_abundances.size(); ++i) {
      (*of) << group_names[i] << '\t' << this->relative_abundances[i] << '\n';
    }
    of->flush();
  } else {
    throw std::runtime_error("Can't write to abundances file.");
  }
}

}

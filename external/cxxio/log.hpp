/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#ifndef CXXIO_LOG_HPP
#define CXXIO_LOG_HPP

#include <chrono>
#include <fstream>
#include <exception>
#include <string>

#include "file.hpp"

class Log : public File::Out {
 public:
  bool verbose;
  std::chrono::time_point<std::chrono::system_clock> start_time;

  Log(std::ostream &stream, bool verbose = true) : File::Out(stream), verbose(verbose) {
    start_time = std::chrono::system_clock::now();
  }
  void flush() {
    if (verbose) {
      std::chrono::time_point<std::chrono::system_clock> end_time = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end_time - start_time;
      std::time_t end = std::chrono::system_clock::to_time_t(end_time);
      std::string end_s(std::ctime(&end));
      end_s.pop_back();
      *this << end_s << " elapsed_time: " << elapsed_seconds.count() << "s\n";
    }
  }
};

template <typename T>
Log& operator<<(Log &os, T t) {
  if (os.verbose) {
    std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
    std::time_t now_t = std::chrono::system_clock::to_time_t(now);
    std::string now_s(std::ctime(&now_t));
    now_s.pop_back();
    os.stream() << now_s << ' ' << t;
    if (!os.stream().good() && !os.stream().eof()) {
      throw std::runtime_error("Error writing type: " + std::string(typeid(T).name()) + " to file " + os.filename());
    }
  }
  return os;
}

#endif

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#ifndef CXXIO_LOG_HPP
#define CXXIO_LOG_HPP

#include <chrono>
#include <fstream>
#include <exception>
#include <string>

#include "cxxio.hpp"
#include "mpi_config.hpp"

class Log : public cxxio::Out {
 public:
  bool verbose;
  bool log_time;
  std::chrono::time_point<std::chrono::system_clock> start_time;

  Log(std::ostream &_stream, bool _verbose = true, bool _log_time = true) : cxxio::Out(_stream), verbose(_verbose), log_time(_log_time) {
    if (this->log_time) {
      start_time = std::chrono::system_clock::now();
    }
  }
  void flush() {
    if (verbose && log_time) {
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
  int rank = 0;
#if defined(MSWEEP_MPI_SUPPORT) && (MSWEEP_MPI_SUPPORT) == 1
  // Only root will log
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
  if (os.verbose && rank == 0) {
    if (os.log_time){
      std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
      std::time_t now_t = std::chrono::system_clock::to_time_t(now);
      std::string now_s(std::ctime(&now_t));
      now_s.pop_back();
      os.stream() << now_s << ' ';
    }
    os.stream () << t;
    if (!os.stream().good() && !os.stream().eof()) {
      throw std::runtime_error("Error writing type: " + std::string(typeid(T).name()) + " to file " + os.filename());
    }
  }
  return os;
}

#endif

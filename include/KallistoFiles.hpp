#ifndef MSWEEP_KALLISTO_FILES_H
#define MSWEEP_KALLISTO_FILES_H

#include <sys/stat.h>

#include <memory>
#include <fstream>
#include <string>

#include "zstr.hpp"

class KallistoFiles {
 private:
  bool file_exists (const std::string& name) const {
    struct stat buffer;   
    return (stat (name.c_str(), &buffer) == 0); 
  }

  std::unique_ptr<zstr::ifstream> open_file(const std::string &path) const {
    if (!file_exists(path)) {
      throw std::runtime_error("File: " + path + " does not exist.");
    }
    zstr::ifstream stream(path);
    if (!stream.good()) {
      throw std::runtime_error("Cannot read from file: " + path + ".");
    }
    return std::unique_ptr<zstr::ifstream>(&stream);
  }

 public:
  KallistoFiles() = default;
  std::unique_ptr<std::istream> ec;
  std::unique_ptr<std::istream> tsv;
  std::unique_ptr<std::istream> cells;
  std::unique_ptr<std::istream> run_info;
  bool batch_mode = false;

  KallistoFiles(std::string path, bool batch_mode) : batch_mode(batch_mode) {
    std::string alignment_path = path + (batch_mode ? "/pseudoalignments" : "/matrix");
    if (batch_mode) {
      cells = open_file(path + "/matrix.cells");
    }
    ec = open_file(alignment_path + ".ec");
    tsv = open_file(alignment_path + ".tsv");
    run_info = open_file(path + "/run_info.json");
  }
};

#endif

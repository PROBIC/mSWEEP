#ifndef MSWEEP_KALLISTO_FILES_HPP
#define MSWEEP_KALLISTO_FILES_HPP

#include <memory>
#include <fstream>
#include <string>

#include "bxzstr.hpp"
#include "cxxio.hpp"

class KallistoFiles {
 private:
  bool file_exists (const std::string& name) const {
    cxxio::In test(name);
    return (!test.stream().fail()); 
  }

  void open_file(const std::string &path, std::unique_ptr<std::istream> &ptr) const {
    if (!file_exists(path)) {
      throw std::runtime_error("File: " + path + " does not exist.");
    }
    ptr.reset(new bxz::ifstream(path));
    if (!ptr->good()) {
      throw std::runtime_error("Cannot read from file: " + path + ".");
    }
  }

 public:
  KallistoFiles() = default;
  std::unique_ptr<std::istream> ec;
  std::unique_ptr<std::istream> tsv;
  std::unique_ptr<std::istream> cells;
  std::unique_ptr<std::istream> run_info;
  bool batch_mode = false;

  KallistoFiles(std::string path, bool batch_mode) : batch_mode(batch_mode) {
    std::string alignment_path = path + (batch_mode ? "/matrix" : "/pseudoalignments");
    if (batch_mode) {
      open_file(path + "/matrix.cells", this->cells);
    }
    open_file(alignment_path + ".ec", this->ec);
    open_file(alignment_path + ".tsv", this->tsv);
    open_file(path + "/run_info.json", this->run_info);
  }
};

#endif

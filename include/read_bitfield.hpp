#ifndef READ_BITFIELD_H
#define READ_BITFIELD_H

#include <memory>
#include <string>
#include <unordered_map>
#include <fstream>

#include "matrix.hpp"
#include "Sample.hpp"
#include "Reference.hpp"
#include "zstr.hpp"

struct KallistoFiles {
  std::unique_ptr<std::istream> ec;
  std::unique_ptr<std::istream> tsv;
  std::unique_ptr<std::istream> cells;
  std::unique_ptr<std::istream> run_info;
  bool batch_mode = false;

  KallistoFiles(std::string path, bool batch_mode) : batch_mode(batch_mode) {
    std::string alignment_path = path + (batch_mode ? "/pseudoalignments" : "/matrix");
    if (batch_mode) {
      cells = std::unique_ptr<std::istream>(new zstr::ifstream(path + "/matrix.cells"));
    }
    ec = std::unique_ptr<std::istream>(new zstr::ifstream(alignment_path + ".ec"));
    tsv = std::unique_ptr<std::istream>(new zstr::ifstream(alignment_path + ".tsv"));
    run_info = std::unique_ptr<std::istream>(new zstr::ifstream(path + "/run_info.json"));
  }
};

void ReadClusterIndicators(std::istream &indicators_file, Reference &reference);
void ReadBitfield(KallistoFiles &kallisto_files, unsigned n_refs, std::vector<Sample> &batch);
void VerifyGrouping(std::istream &run_info, unsigned n_refs);

#endif

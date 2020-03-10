#ifndef FILE_H
#define FILE_H

#include <string>
#include <fstream>
#include <memory>
#include <exception>
#include <iostream>

#include "bxzstr.hpp"

namespace File {
  namespace exceptions {
    struct file_exception : public std::exception {
      std::string msg;
      virtual ~file_exception() override = default;
      file_exception(std::string msg) : msg(msg) {}
      const char* what() const noexcept override { return msg.c_str(); }
    };
    struct file_not_writable : public file_exception {
      file_not_writable(std::string name) : file_exception("File " + name + " is not writable (does the directory exist?).") {}
    };
    struct cannot_read_from_file : public file_exception {
      cannot_read_from_file(std::string name) : file_exception("Cannot read from file: " + name + ".") {}
    };
  }

  class Out {
    std::unique_ptr<std::ostream> byname;
    std::ostream os;
    std::string myfile;
  public:
    Out(std::ostream& stream = std::cout) : byname(), os(stream.rdbuf()) { }
    Out(std::string filename) : byname(new std::ofstream(filename)), os(this->byname->rdbuf()), myfile(filename) {
      if (!os)
	throw exceptions::file_not_writable(filename);
    }
    void open(const std::string &filename) {
      os.flush();
      byname.reset(new std::ofstream(filename));
      os.rdbuf(byname->rdbuf());
      os.setstate(byname->rdstate());
      if (!os)
	throw exceptions::file_not_writable(filename);
    }
    void open_compressed(const std::string &filename) {
      os.flush();
      byname.reset(new bxz::ofstream(filename));
      os.rdbuf(byname->rdbuf());
      os.setstate(byname->rdstate());
      if (!os)
	throw exceptions::file_not_writable(filename);
    }      
    void close() {
      os.flush();
      byname.reset();
      os.rdbuf(std::cout.rdbuf());
      os.setstate(std::cout.rdstate());
    }
    std::ostream& stream() { return this->os; }
    const std::string& filename() const { return this->myfile; }
  };

  class In {
    std::unique_ptr<bxz::ifstream> byname;
    std::istream is;
    std::string myfile;

    void reset_state(const std::string &filename) {
      byname.reset(new bxz::ifstream(filename));
      is.rdbuf(byname->rdbuf());
      is.setstate(byname->rdstate());      
    }
  public:
    In(std::istream& stream = std::cin) : byname(), is(stream.rdbuf()) { }
    In(std::string filename) : is(std::cout.rdbuf()), myfile(filename) {
      this->open(filename);
    }
    In(In&& in) : is(std::cout.rdbuf()) {
      if (!in.filename().empty()) {
	this->open(in.filename());
      }
    }
    void open(const std::string &filename) {
      myfile = filename;
      reset_state(myfile);
      if (!is)
      	throw std::runtime_error("Cannot read from file: " + filename + ".");
    }
    void close() {
      byname.reset();
    }
    void rewind() {
      reset_state(myfile);
    }
    std::istream& stream() { return this->is; }
    const std::string& filename() const { return this->myfile; }
  };
  template <typename T>
  Out& operator<<(Out &os, T t) {
    os.stream() << t;
    if (!os.stream().good() && !os.stream().eof()) {
      throw std::runtime_error("Error writing type: " + std::string(typeid(T).name()) + " to file " + os.filename());
    }
    return os;
  }
  template <typename T>
  In& operator>>(In &is, T &t) {
    is.stream() >> t;
    if (!is.stream().good() && !is.stream().eof()) {
      throw std::runtime_error("Error writing type: " + std::string(typeid(T).name()) + " to file " + is.filename());
    }
    return is;
  }
}

#endif

// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 26/05/2005

#ifndef NameReader_HPP
#define NameReader_HPP
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <iostream>
#include <CoinHelperFunctions.hpp>
#include "OsiSolverInterface.hpp"
//#include <tr1/unordered_map>
#include <map>

namespace Bonmin{
/** A class for reading a .col or .row file containing name for variables and constraints (usually ampl generated file).
   */
class NamesReader
{
public:
  /** Constructor with a file name given by a const char * */
  NamesReader(const char * fileName, const char * suffix);
  /** Constructor with a file name given by a string and also default (empty string) */
  NamesReader(const std::string & fileName="", const std::string& suffix=".col");
  /** Reads the .col file*/
  bool readFile();
  /** Reads the .col file fileName*/
  bool readFile(const std::string &file)
  {
    file_=file;
    return readFile();
  }

  /** Copy the names to Names. */
  void copyNames(OsiSolverInterface::OsiNameVec& Names);

  /** Access Names of indexed by i. */
  const std::string& name(int i){
   return names_[i];
  }

  /** Access index of variable str */
  int index(const char * str){
    return indices_[str];
  }
private:
  /// Name of the file to read.
  std::string file_;

  /// Suffix of the file (".col", ".row")
  std::string suffix_;

  /// String comparison strucutre.
  struct ltstr
  {
    bool operator()(const char* s1, const char* s2) const
    {
      return strcmp(s1, s2) < 0;
    }
  };

  /// Hash type.
  typedef std::map<const char *, int, ltstr> namesHash;
  ///Hash map used to store the indices.
  namesHash indices_;
  ///Variable names.
  std::vector<std::string> names_;
};
}
#endif

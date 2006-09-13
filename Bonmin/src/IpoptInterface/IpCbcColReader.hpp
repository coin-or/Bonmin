// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 26/05/2005

#ifndef IpCbcColReader_HPP
#define IpCbcColReader_HPP
#include <string>
#include <list>
#include <fstream>
#include <iostream>
#include <CoinHelperFunctions.hpp>
#include <map>
#include <vector>

/** A class for reading a .col file containing name for variable (ampl generated file).
   */
class IpCbcColReader
{
public:
  /** Constructor with a file name given by a const char * */
  IpCbcColReader(const char * fileName);
  /** Constructor with a file name given by a string and also default (empty string) */
  IpCbcColReader(const std::string & fileName="");
  /** Reads the .col file*/
  bool readFile();
  /** Reads the .col file fileName*/
  bool readFile(const std::string &fileName)
  {
    fileName_=fileName;
    return readFile();
  }

  /** Copy the names to varNames */
  void copyNames(std::string *varNames, int n_var);
private:
  /// Name of the file to read.
  std::string fileName_;

  /// String comparison strucutre.
  struct ltstr
  {
    bool operator()(const char* s1, const char* s2) const
    {
      return strcmp(s1, s2) < 0;
    }
  };

  /// Hash type.
  //typedef __gnu_cxx::hash_map<const char *, int> namesHash;//, __gnu_cxx::hash <const char *>, eqstr > namesHash;
  typedef std::map< const char *, int, ltstr > namesStorage;

  ///Hash map used to store the indices.
  namesStorage varIndices_;
  ///Variable names.
  std::vector<std::string> varNames_;
};

#endif

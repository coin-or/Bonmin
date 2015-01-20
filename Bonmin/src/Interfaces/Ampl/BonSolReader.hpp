// (C) Copyright CNRS 2011
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, LIF, CNRS,
//
// Date : 03/01/2011

#ifndef SolReader_HPP
#define SolReader_HPP
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <iostream>
#include <CoinHelperFunctions.hpp>
#include "OsiSolverInterface.hpp"
#include "BonTypes.hpp"

namespace Bonmin{
/** A class for reading a .col or .row file containing name for variables and constraints (usually ampl generated file).
   */
class SolReader
{
public:
  /** Constructor with a file name given by a const char * */
  SolReader(const char * fileName, const char * suffix);
  /** Constructor with a file name given by a string and also default (empty string) */
  SolReader(const std::string & fileName="", const std::string& suffix=".col");
  /** Reads the .sol file*/
  bool readFile();
  /** Reads the .sol file fileName*/
  bool readFile(const std::string &file)
  {
    file_=file;
    return readFile();
  }

  /** Copy the names to Names. */
  void copySol(double * x);

  const double * x(){
    return x_();
  }

  /** Set the number of variables in the problem.*/
  void set_n_cols(int n){
    x_.resize(n);
  }
private:
  /// Name of the file to read.
  std::string file_;

  /// Suffix of the file (".col", ".row")
  std::string suffix_;

  /// Sol values
 vector<double> x_;
};
}
#endif

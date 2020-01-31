// (C) Copyright CNRS 2011
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, LIF, CNRS,
//
// Date : 03/01/2011
#include "BonSolReader.hpp"
#include <fstream>
#include <iostream>

namespace Bonmin{

SolReader::SolReader(const char * file, const char * suffix)
    :
    file_(), suffix_(suffix), x_()
{
  assert(file!= NULL);
  file_=file;
  if (suffix!=NULL)
    suffix_ = suffix;
}
SolReader::SolReader(const std::string & file, const std::string & suffix)
    :
    file_(file), suffix_(suffix), x_()
{}

bool SolReader::readFile()
{
  std::string fileName = file_;
  size_t size = fileName.size();
  bool hasNlExtension =  (fileName.size()>4) && (fileName[size - 1] =='l') && (fileName[size - 2] =='n') && (fileName[size - 3] =='.');
  if(hasNlExtension)
    fileName.erase(size-3,3);
  fileName+=suffix_;
  std::ifstream inFile(fileName.c_str());
  if(!inFile.is_open()) {
    return false;
  }
  std::string token;
  inFile>>token;
  assert(token == "bonmin:");

  std::string status;
  inFile>>status;
  inFile>>token;
  if(token == "Options"){
    for(int i = 0 ; i < 6 ; i++){
      inFile>>token;
    }
    int n_cols, n_cols_2;
    inFile>>n_cols_2>>n_cols;
    if(n_cols != static_cast<int>(x_.size())){
       fprintf(stderr, "Number of columns different %d\n", n_cols);
       x_.resize(n_cols);
    }
  }
for(size_t i = 0 ; i < x_.size() ; i++){
     inFile>>x_[i];
  }
  return true;
}

void
SolReader::copySol(double * x)
{
  std::copy(x_.begin(), x_.end(), x);
}
}

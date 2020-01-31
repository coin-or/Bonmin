// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 26/05/2005
#include "BonColReader.hpp"
#include <fstream>
#include <iostream>

namespace Bonmin{

NamesReader::NamesReader(const char * file, const char * suffix)
    :
    file_(), suffix_(suffix), indices_(), names_()
{
  assert(file!= NULL);
  file_=file;
  if (suffix!=NULL)
    suffix_ = suffix;
}
NamesReader::NamesReader(const std::string & file, const std::string & suffix)
    :
    file_(file), suffix_(suffix), indices_(), names_()
{}

bool NamesReader::readFile()
{
  std::string colFileName = file_;
  size_t size = colFileName.size();
  bool hasNlExtension =  (colFileName.size()>4) && (colFileName[size - 1] =='l') && (colFileName[size - 2] =='n') && (colFileName[size - 3] =='.');
  if(hasNlExtension)
    colFileName.erase(size-3,3);
  colFileName+=suffix_;
  std::ifstream inFile(colFileName.c_str());
  if(!inFile.is_open()) {
    return false;
  }
  std::string name;
  int nVar = 0;
  do {
    name="";
    inFile>>name;
    if(name.size()==0)
      continue;
    names_.push_back(name);
    indices_[names_[nVar].c_str()] = nVar;
    nVar++;
  }
  while(!inFile.eof());

  //  names_ = new std::string[nVar];
  for(int i = 0 ; i < nVar ; i++) {
    assert(i==indices_ [ names_ [i].c_str()]);
  }
  return true;
}

void
NamesReader::copyNames(OsiSolverInterface::OsiNameVec& names)
{
  names_ = names;
}
}

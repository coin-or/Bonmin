// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 26/05/2005

#include "BonBoundsReader.hpp"

#include <fstream>

namespace Bonmin {
BoundsReader::~BoundsReader()
{
  gutsOfDestructor();
}
void BoundsReader::gutsOfDestructor()
{
  if(nLower_ > 0) {
    assert(lowerBounds_!= NULL);
    delete [] lowerBounds_;
    lowerBounds_ = NULL;
    assert(indexLowers_ != NULL);
    delete [] indexLowers_;
    indexLowers_ = NULL;
  }
  else {
    assert(lowerBounds_ == NULL);
    assert(indexLowers_ == NULL);
  }
  if(nUpper_ > 0) {
    assert(upperBounds_!= NULL);
    delete [] upperBounds_;
    upperBounds_ = NULL;
    assert(indexUppers_ != NULL);
    delete [] indexUppers_;
    indexUppers_ = NULL;
  }
  else {
    assert(upperBounds_ == NULL);
    assert(indexUppers_ == NULL);
  }
  nLower_=0;
  nUpper_=0;
}

void BoundsReader::read(const std::string &fileName)
{
  setFileName(fileName);
  read();
}

void BoundsReader::read()
{
  gutsOfDestructor();
  //First count the number of lower and upper bounds resp
  std::string lo="LO";
  std::string up="UP";
  std::string in;
  {
    std::ifstream fin(fileName_.c_str());
    //std::streampos begin = fin.tellg();
    while(!fin.eof()) {
      fin>>in;
      if(in==lo)
        nLower_++;
      else if(in==up)
        nUpper_++;
      else
        throw;
      fin.ignore(10000,'\n');
    }
  }
  if(nLower_ > 0) {
    lowerBounds_ = new double[nLower_];
    indexLowers_ = new int[nLower_];
  }
  if(nUpper_ > 0) {
    upperBounds_ = new double[nUpper_];
    indexUppers_ = new int[nUpper_];
  }
  //fin.seekg(0);
  // fin.seekg( 0, std::ios::beg);
  nLower_ = 0;
  nUpper_ = 0;
  {
    std::ifstream fin2(fileName_.c_str());
    //fin.close();
    //fin.open(fileName_.c_str());
    while(!fin2.eof()) {
      int index;
      double bound;
      fin2>>in>>index>>bound;
      if(in==lo) {
        lowerBounds_[nLower_] = bound;
        indexLowers_[nLower_++] = index;
      }
      else if(in==up) {
        upperBounds_[nUpper_] = bound;
        indexUppers_[nUpper_++] = index;
      }
      else
        throw;
      fin2.ignore(10000,'\n');
    }
  }
}

void BoundsReader::readAndApply(OsiTMINLPInterface * solver)
{
  read();
  for(int i = 0 ; i < nLower_ ; i++) {
    solver->setColLower(indexLowers_[i], lowerBounds_[i]);
  }
  for(int i = 0 ; i < nUpper_ ; i++) {
    solver->setColUpper(indexUppers_[i], upperBounds_[i]);
  }
}
}

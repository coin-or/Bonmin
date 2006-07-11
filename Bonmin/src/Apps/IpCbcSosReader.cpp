// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 05/26/2005

#include "IpCbcSosReader.hpp"
#include <fstream>
#include <iostream>
/** Keywords of the .int file */
static std::string keywords[]= { "priorities", "PRIORITIES",
                                 "special_ordered_set_1",
                                 "SPECIAL_ORDERED_SET_1"};
IpCbcSosReader::IpCbcSosReader(const char * fileName)
    :
    fileName_(), sos_()
{
  assert(fileName != NULL);
  fileName_=fileName;
}
IpCbcSosReader::IpCbcSosReader(const std::string & fileName)
    :
    fileName_(fileName), sos_()
{}

void IpCbcSosReader::readFile()
{
  std::ifstream inFile(fileName_.c_str());
  if(!inFile.is_open()) {
    std::cerr<<"Error in opening SOS file";
  }
  std::string keyword;
  do {
    keyword="";
    inFile>>keyword;
    if(keyword.size()==0)
      continue;
    inFile.ignore(1024,'\n');
    if(keyword == keywords[0] || keyword ==keywords[1])
      readPrioritiesSection(inFile);
    else if(keyword == keywords[2] || keyword ==keywords[3])
      readSosSection(inFile);
    else
      std::cerr<<"unrecognized keyword "<<keyword;
  }
  while(!inFile.eof());
  print();
}

void IpCbcSosReader::readPrioritiesSection(std::ifstream & inFile)
{
  //For the moment skip priorities on variable
  std::string key;
  do {
    inFile>>key;
    if(key=="END")
      break;
    else
      inFile.ignore(1024,'\n');
  }
  while(1);
  inFile>>key;
  assert(key==keywords[0] || key == keywords[1]);
}

void IpCbcSosReader::readSosSection(std::ifstream & inFile)
{
  std::string key;
  do {
    inFile>>key;
    if(key!="END") {
      int sosPriority;
      int sosNElements;
      inFile>>sosPriority>>sosNElements;
      inFile.ignore(1024,'\n');
      sos_.push_back(sosDesc(key,sosPriority,sosNElements));
      for(int i=0 ; i < sosNElements ; i++) {
        inFile>>sos_.back().var_[i]>>sos_.back().varPriority_[i];
        inFile.ignore(1024,'\n');
      }
    }
    else
      break;
  }
  while(1);
  inFile>>key;
  inFile.ignore(1024,'\n');
  assert(key==keywords[2] || key == keywords[3]);
}
void IpCbcSosReader::print()
{
  for(std::list<sosDesc>::iterator i = sos_.begin(); i!=sos_.end(); i++)
    i->print();
}

void IpCbcSosReader::readNames()
{
  std::string colFileName = fileName_;
  int size = colFileName.size();
  assert(colFileName[size - 1] =='t');
  assert(colFileName[size - 2] =='n');
  assert(colFileName[size - 3] =='i');
  assert(colFileName[size - 4] =='.');
  colFileName.erase(size-4,4);
  colFileName+=".col";
  std::ifstream inFile(colFileName.c_str());
  if(!inFile.is_open()) {
    std::cerr<<"Error in opening Names file";
    return;
  }
  std::string name;
  int nVar = 0;
  do {
    name="";
    inFile>>name;
    if(name.size()==0)
      continue;
    varNames_[name.c_str()]=nVar++;
  }
  while(!inFile.eof());
}

/** Add the sos constraints to the model */
void addSosToModel(CbcModel & model)
{}

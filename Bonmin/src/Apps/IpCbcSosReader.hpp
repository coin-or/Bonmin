// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 05/26/2005

#ifndef IpCbcSosReader_HPP
#define IpCbcSosReader_HPP
#include <string>
#include <list>
#include <fstream>
#include <iostream>
#include <CoinHelperFunctions.hpp>
#include <ext/hash_map>
#include <CbcModel.hpp>
/** A class for reading a .int file containing sos constraints in Leyffer's
 *  format
 *  */
class IpCbcSosReader
{
public:
  /** Constructor with a file name given by a const char * */
  IpCbcSosReader(const char * fileName);
  /** Constructor with a file name given by a string and also default (empty string) */
  IpCbcSosReader(const std::string & fileName="");
  /** Reads the .int file*/
  void readFile();
  /** Reads the .int file fileName*/
  void readFile(const std::string &fileName)
  {
    fileName_=fileName;
    readFile();
  }
  /** Read the priority section of a .int file */
  void readPrioritiesSection(std::ifstream & inFile);
  /** Read the SOS section section of a .int file */
  void readSosSection(std::ifstream & inFile);
  /** Read variable names in a .col file */
  void readNames();
  /** Add the sos constraints to the model */
  void addSosToModel(CbcModel & model);
  /** print the SOS constraints */
  void print();
private:
  struct sosDesc
  {
    std::string name_;
    int nElem_;
    double sosPriority_;
    std::string * var_;
    double * varPriority_;
    sosDesc(const std::string &name, double sosPriority, int nElem):
        name_(name), sosPriority_(sosPriority), nElem_(nElem),
        var_(NULL), varPriority_(NULL)
    {
      var_ = new std::string [ nElem_ ];
      varPriority_ = new double [ nElem_ ];
    }
    sosDesc(const sosDesc &toCopy):
        name_(toCopy.name_), sosPriority_(toCopy.sosPriority_),
        nElem_(toCopy.nElem_),
        var_(NULL), varPriority_(NULL)
    {
      var_ = new std::string [ nElem_ ];
      CoinCopyN(toCopy.var_, nElem_, var_);
      varPriority_ = new double [ nElem_ ];
      CoinCopyN(toCopy.varPriority_, nElem_, varPriority_);
    }
    ~sosDesc()
    {
      delete [] var_;
      delete [] varPriority_;
    }
    void print()
    {
      std::cout<<name_<<" "<<nElem_<<" "<<sosPriority_<<std::endl;
      for(int i=0 ; i < nElem_ ; i++) {
        std::cout<<var_[i]<<" "<<varPriority_[i]<<std::endl;
      }
    }
private :
    sosDesc()
    {}
  }
  ;
  std::string fileName_;
  std::list<sosDesc> sos_;
  struct eqstr
  {
    bool operator()(const char* s1, const char* s2) const
    {
      return strcmp(s1, s2) == 0;
    }
  };
  __gnu_cxx::hash_map<const char *, int, __gnu_cxx::hash <const char *>, eqstr > varNames_;
};

#endif

// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 26/05/2005

#ifndef _IPCBCSTARTPOINTREADER_H_
#define _IPCBCSTARTPOINTREADER_H_
#include <string>
#include <list>
#include <fstream>
#include <iostream>
#include "IpoptInterface.hpp"

/** This class reads a file with a starting point for Ipopt initalization. File format is number of primals number of duals then values one after another
 * Numbering of variables is first variables, then duals on lower bounds duals on upper bounds and to finish duals on constraints */
class IpCbcStartPointReader
{
public:
  /** Constructor with fileName_ given by a string (and default) */
  IpCbcStartPointReader(std::string fileName = ""):
      fileName_(fileName),
      primals_(NULL),
      duals_(NULL)
  {}
  /** Constructor with fileName_ given by a const char * */
  IpCbcStartPointReader(const char * fileName):
      fileName_(fileName),
      primals_(NULL),
      duals_(NULL)
  {}

  /** Reads the .initP file*/
  bool readFile();
  /** Reads the .initP file fileName*/
  bool readFile(const std::string &fileName)
  {
    fileName_=fileName;
    return readFile();
  }
  /** Read warmstart info and apply to an IpoptInterface */
  bool readAndApply(IpoptInterface& solver);
  ~IpCbcStartPointReader()
  {
    gutsOfDestructor();
  }


  /// Dealocate arrays
  void gutsOfDestructor()
  {
    if(primals_!=NULL)
      delete[] primals_;
    if(duals_!=NULL)
      delete[] duals_;
  }

  /// Access primal variables values.
  const double * getPrimals()
  {
    return primals_;
  }
  /// Access dual variables values.
  const double * getDuals()
  {
    return duals_;
  }
private:
  /** Name of the file with initial point */
  std::string fileName_;

  /// Primal variables values.
  double * primals_;
  /// Dual variables values.
  double * duals_;
};
#endif /*_IPCBCINITPOINTREADER_H_*/

// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 26/05/2005

#ifndef BoundsReader_HPP
#define BoundsReader_HPP

#include <string>
#include "BonOsiTMINLPInterface.hpp"


namespace Bonmin {
/** Reads a file containing change bounds for variables.
    Files follows pretty much the Bounds section in MPS standard.*/
class BoundsReader
{
public:
  //Default constructor
  BoundsReader():
      fileName_(),
      lowerBounds_(NULL),
      upperBounds_(NULL),
      nLower_(0),
      nUpper_(0)
  {}

  // Constructor with name of the file to read passed.
  BoundsReader(const std::string &fileName):
      fileName_(fileName),
      lowerBounds_(NULL),
      upperBounds_(NULL),
      indexLowers_(NULL),
      indexUppers_(NULL),
      nLower_(0),
      nUpper_(0)
  {}

  // Set the name of the file to read.
  void setFileName(const std::string &fileName)
  {
    fileName_ = fileName;
  }

  // Destructor
  ~BoundsReader();

  // Cleanup allocated data
  void gutsOfDestructor();


  // Read the file with given fileName
  void read(const std::string &);

  //Read the file named fileName_
  void read();

  //Read fileName_ and apply the bounds read to solver
  void readAndApply(OsiTMINLPInterface * solver);
private:

  /// Current file
  std::string fileName_;
  /// changed lower bounds
  double * lowerBounds_;
  /// changed upper bounds
  double * upperBounds_;
  /// index of the changed lowerBounds_
  int * indexLowers_;
  /// index of the changed upperBounds_
  int * indexUppers_;
  /// number of changed lowerBounds_
  int nLower_;
  /// number of changed upperBounds_
  int nUpper_;
};
}
#endif

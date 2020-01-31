// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 26/05/2005
#include "BonStartPointReader.hpp"


namespace Bonmin {

  bool StartPointReader::readFile()
  {
    std::ifstream inFile(fileName_.c_str());
    if(!inFile.is_open()) {
      std::cerr<<"Error in opening initial point file";
      return false;
    }
    int numPrimals;
    int numDuals;
    inFile>>numPrimals>>numDuals;
    gutsOfDestructor();
    primals_ = new double [numPrimals];
    duals_ = new double[numDuals];
    for(int i = 0; i < numPrimals ; i++) {
      inFile>>primals_[i];
    }
    for(int i = 0; i < numDuals ; i++) {
      inFile>>duals_[i];
    }
    return true;
  }

  bool StartPointReader::readAndApply(OsiTMINLPInterface * solver)
  {
    readFile();
    solver->solver()->enableWarmStart();
    if(primals_)
      solver->setColSolution(primals_);
    else {
      std::cerr<<"No warm start info ???"<<std::endl;
      return 0;
    }
    if(duals_)
      solver->setRowPrice(duals_);
    else {
      std::cerr<<"No warm start info ???"<<std::endl;
      return 0;
    }
    return 1;
  }
}

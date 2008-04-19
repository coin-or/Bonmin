// (C) Copyright CNRS 2008
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, LIF, CNRS-Marseille Universites.
//
// Date : 03/17/2008
#ifndef BonRandomChoice_H
#define BonRandomChoice_H

#include "OsiChooseVariable.hpp"
#include "OsiSolverInterface.hpp"

class BonRandomChoice : public OsiChooseVariable {
  public:
  ///Default constructor
  BonRandomChoice(): OsiChooseVariable(){
  }

  //Constructor from solver
  BonRandomChoice(const OsiSolverInterface * solver):
    OsiChooseVariable(solver){
  }

  // Copy constructor
  BonRandomChoice(const BonRandomChoice &other):
    OsiChooseVariable(other){
  }

  // Assignment operator
  BonRandomChoice & operator=(const BonRandomChoice &rhs){
    OsiChooseVariable::operator=(rhs);
    return (*this);
  }

  // Virtual copy
  virtual OsiChooseVariable * clone() const{
     return new BonRandomChoice(*this);
  }

  /// Destructor
  virtual ~BonRandomChoice(){
  }

  virtual int chooseVariable( OsiSolverInterface * solver,
                              OsiBranchingInformation * info,
                              bool fixVariables){
    if(numberUnsatisfied_){
      int chosen = (int) (floor(CoinDrand48() * (numberUnsatisfied_)));
      bestObjectIndex_ = list_[chosen];
      bestWhichWay_ = solver->object(bestObjectIndex_)->whichWay();
      firstForcedObjectIndex_ = -1;
      firstForcedWhichWay_ =-1;
      return 0;
    }
    else {
      return 1;
    }
  }
};
#endif

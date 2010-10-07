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
#include <list>

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

/** Own version of setupList since Osi version is broken and what we want to do here is anyway much simpler.*/

  virtual int setupList(OsiBranchingInformation * info, bool initialize){
  if (initialize) {
    status_=-2;
    delete [] goodSolution_;
    bestObjectIndex_=-1;
    numberStrongDone_=0;
    numberStrongIterations_ = 0;
    numberStrongFixed_ = 0;
    goodSolution_ = NULL;
    goodObjectiveValue_ = COIN_DBL_MAX;
    numberOnList_=0;
    }
    numberUnsatisfied_=0;
    int numberObjects = solver_->numberObjects();
    assert (numberObjects);
    int bestPriority=COIN_INT_MAX;
    std::fill(list_, list_+numberObjects, -1);
    OsiObject ** object = info->solver_->objects();
    // Say feasible
    bool feasible = true;
    for (int i=0;i<numberObjects;i++) {
    int way;
    double value = object[i]->infeasibility(info,way);
    if (value>0.0) {
      numberUnsatisfied_++;
      if (value==COIN_DBL_MAX) {
	// infeasible
	feasible=false;
	break;
      }
      int priorityLevel = object[i]->priority();
      // Better priority? Flush choices.
      if (priorityLevel<bestPriority) {
        numberUnsatisfied_ = 0;
        std::fill(list_, list_+numberObjects, -1);
      }
      bestPriority = priorityLevel;
      if (priorityLevel==bestPriority) {
        list_[numberUnsatisfied_]=i;
        numberUnsatisfied_++;
        }
      }
    }
    return numberUnsatisfied_;
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

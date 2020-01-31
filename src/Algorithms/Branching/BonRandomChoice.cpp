// (C) Copyright CNRS 2008
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, LIF, CNRS-Marseille Universites.
//
// Date : 03/17/2008
#include "BonRandomChoice.hpp"
#include "BonminConfig.h"
#include "CoinFinite.hpp"

int BonRandomChoice::setupList(OsiBranchingInformation * info, bool initialize){
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
    //bool feasible = true;
    for (int i=0;i<numberObjects;i++) {
    int way;
    double value = object[i]->infeasibility(info,way);
    if (value>0.0) {
      if (value==COIN_DBL_MAX) {
  // infeasible
  //feasible=false;
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

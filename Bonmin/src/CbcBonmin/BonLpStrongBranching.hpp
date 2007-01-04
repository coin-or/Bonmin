// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef BonLpStrongBranching_H
#define BonLpStrongBranching_H

#include "BonChooseVariable.hpp"

namespace Bonmin {
  class LpStrongBranching : public BonChooseVariable {

 public:
  

 /// Constructor from a solver
 LpStrongBranching(OsiTMINLPInterface * solver);

 /// Copy constructor
 LpStrongBranching(const LpStrongBranching &other);

 ///Assignment operator
 LpStrongBranching & operator=(const LpStrongBranching &rhs);

 ///Virtual copy constructor
 virtual OsiChooseVariable * clone() const;

 /// Destructor
 virtual ~LpStrongBranching();

 void setMaxCuttingPlaneIter_(int num){
   maxCuttingPlaneIteration_ = num;
 }

  /// This method determines the predicted changes in the objective
  /// function. It has to be implemented by any class that inherits
  /// off of this class.
  virtual int fill_changes(OsiSolverInterface * solver,
			   OsiBranchingInformation *info,
			   bool fixVariables,
			   int numStrong, double* change_down,
			   double* change_up, int& best_way);
 private:
 int maxCuttingPlaneIteration_;
};

}
#endif

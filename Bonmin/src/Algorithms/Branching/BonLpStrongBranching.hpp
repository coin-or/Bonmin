// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef BonLpStrongBranching_H
#define BonLpStrongBranching_H

#include "BonChooseVariable.hpp"

namespace Bonmin {
  class LpStrongBranching : public BonChooseVariable {

 public:
    /** Default constructor. */
 LpStrongBranching(OsiTMINLPInterface *si);


    //Copy constructor
 LpStrongBranching(const LpStrongBranching &other);

 ///Assignment operator
 LpStrongBranching & operator=(const LpStrongBranching &rhs);

 ///Virtual copy constructor
 virtual OsiChooseVariable * clone() const;

 /// Destructor
 virtual ~LpStrongBranching();

 void setMaxCuttingPlaneIter(int num){
   maxCuttingPlaneIterations_ = num;
 }

  /// This method determines the predicted changes in the objective
  /// function. It has to be implemented by any class that inherits
  /// off of this class.
  virtual int fill_changes(OsiSolverInterface * solver,
			   OsiBranchingInformation *info,
			   bool fixVariables,
			   int numStrong, double* change_down,
			   double* change_up, int& best_way);
 
 static void registerOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions);
 private:
 int maxCuttingPlaneIterations_;
};

}
#endif

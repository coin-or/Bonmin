// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef BonQPStrongBranching_H
#define BonQPStrongBranching_H

#include "BonChooseVariable.hpp"
#include "BonBranchingTQP.hpp"

namespace Bonmin {

/** This class chooses a variable to branch on

    This implementation solves the QP model for different branches
    (strong branching).
*/

class BonQPStrongBranching : public BonChooseVariable  {
 
public:

  /// Constructor from solver (so we can set up arrays etc)
  BonQPStrongBranching (OsiTMINLPInterface * solver,
			bool solve_nlp = false);

  /// Copy constructor 
  BonQPStrongBranching (const BonQPStrongBranching &);
   
  /// Assignment operator 
  BonQPStrongBranching & operator= (const BonQPStrongBranching& rhs);

  /// Clone
  virtual OsiChooseVariable * clone() const;

  /// Destructor
  virtual ~BonQPStrongBranching ();

#ifdef BONMIN_CURVATURE_USE_QP_IF_SOS_FOUND
public:
#else
protected:
#endif
  virtual int fill_changes(OsiSolverInterface * solver,
			   OsiBranchingInformation *info,
			   bool fixVariables,
			   int numStrong, double* change_down,
			   double* change_up, int& best_way);  

  bool solve_nlp_;

private:
  /// Default Constructor 
  BonQPStrongBranching ();

};

}

#endif

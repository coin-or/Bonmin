// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef BonCurvBranching_H
#define BonCurvBranching_H

#include "BonChooseVariable.hpp"
#include "BonCurvatureEstimator.hpp"
#include "BonOsiTMINLPInterface.hpp"

namespace Bonmin {

/** Implementation of BonChooseVariable for curvature-based braching.
*/

class BonCurvBranching : public BonChooseVariable  {
 
public:

  /// Constructor from solver (so we can set up arrays etc)
  BonCurvBranching (OsiTMINLPInterface * solver);

  /// Copy constructor 
  BonCurvBranching (const BonCurvBranching &);
   
  /// Assignment operator 
  BonCurvBranching & operator= (const BonCurvBranching& rhs);

  /// Clone
  virtual OsiChooseVariable * clone() const;

  /// Destructor
  virtual ~BonCurvBranching ();

protected:
  virtual int fill_changes(OsiSolverInterface * solver,
			   OsiBranchingInformation *info,
			   bool fixVariables,
			   int numStrong, double* change_down,
			   double* change_up, int& best_way);  

  SmartPtr<CurvatureEstimator> cur_estimator_;

private:
  /// Default Constructor 
  BonCurvBranching ();

};

}

#endif

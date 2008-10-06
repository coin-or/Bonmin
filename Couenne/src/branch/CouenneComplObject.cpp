/*
 * Name:    CouenneComplObject.cpp
 * Authors: Pietro Belotti, Lehigh University
 * Purpose: Implementation of branching rules for complementarity constraints
 *
 * (C) Carnegie-Mellon University, 2008.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneComplObject.hpp"


/// empty constructor (for unused objects)
CouenneComplObject::CouenneComplObject () {

}


/// Constructor with information for branching point selection strategy
CouenneComplObject::CouenneComplObject (CouenneProblem *p, 
					exprVar *ref, Bonmin::BabSetupBase *base, JnlstPtr jnlst):
  CouenneObject (p,  ref, base, jnlst) {

}


/// Constructor with lesser information, used for infeasibility only
CouenneComplObject::CouenneComplObject (exprVar *ref, Bonmin::BabSetupBase *base, JnlstPtr jnlst):
  CouenneObject (ref, base, jnlst) {

}


/// Copy constructor
CouenneComplObject::CouenneComplObject (const CouenneObject &src): 
  CouenneObject (src) {

}
    

/// compute infeasibility of this variable, |w - f(x)| (where w is
/// the auxiliary variable defined as w = f(x)
double CouenneComplObject::infeasibility (const OsiBranchingInformation *info, int &way) const {

  return 0.;
}


/// compute infeasibility of this variable, |w - f(x)|, where w is
/// the auxiliary variable defined as w = f(x)
double CouenneComplObject::checkInfeasibility (const OsiBranchingInformation * info) const {

  return 0.;
}


/// fix (one of the) arguments of reference auxiliary variable 
double CouenneComplObject::feasibleRegion (OsiSolverInterface*, const OsiBranchingInformation*) const {
  return 0;
}


/// create CouenneBranchingObject or CouenneThreeWayBranchObj based
/// on this object
OsiBranchingObject *CouenneComplObject::createBranch (OsiSolverInterface*, 
						      const OsiBranchingInformation*, int) const {
  return NULL;
}

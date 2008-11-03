/*
 * Name:    CouenneComplObject.cpp
 * Authors: Pietro Belotti, Lehigh University
 * Purpose: Implementation of branching rules for complementarity constraints
 *
 * (C) Carnegie-Mellon University, 2008.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneComplObject.hpp"
#include "CouenneComplBranchingObject.hpp"


/// Constructor with information for branching point selection strategy
CouenneComplObject::CouenneComplObject (CouenneProblem *p, 
					exprVar *ref, Bonmin::BabSetupBase *base, JnlstPtr jnlst):
  CouenneObject (p, ref, base, jnlst) {
  jnlst -> Printf (J_DETAILED, J_PROBLEM, "[created Complementarity constraint object]\n");
}


/// Constructor with lesser information, used for infeasibility only
CouenneComplObject::CouenneComplObject (exprVar *ref, Bonmin::BabSetupBase *base, JnlstPtr jnlst):
  CouenneObject (ref, base, jnlst) {}


/// Copy constructor
CouenneComplObject::CouenneComplObject (const CouenneObject &src): 
  CouenneObject (src) {}
    

/// compute infeasibility of this variable, |w - f(x)| (where w is
/// the auxiliary variable defined as w = f(x)
double CouenneComplObject::infeasibility (const OsiBranchingInformation *info, int &way) const {

  expression **arglist = reference_ -> Image () -> ArgList ();

  int index0 = arglist [0] -> Index (),
      index1 = arglist [1] -> Index ();

  CouNumber 
    x0 = fabs (info -> solution_ [index0]),
    x1 = fabs (info -> solution_ [index1]);

  // if x1 < x0, it is preferrable to branch with x1=0 instead of x0=0
  // as this is closer to the point
  way = (x1 < x0) ? 1 : 0;

  return x0 * x1;
}


/// compute infeasibility of this variable, |w - f(x)|, where w is
/// the auxiliary variable defined as w = f(x)
double CouenneComplObject::checkInfeasibility (const OsiBranchingInformation * info) const {

  expression **arglist = reference_ -> Image () -> ArgList ();

  int index0 = arglist [0] -> Index (),
      index1 = arglist [1] -> Index ();

  return 
    fabs (info -> solution_ [index0]) * 
    fabs (info -> solution_ [index1]);
}


/// create CouenneBranchingObject or CouenneThreeWayBranchObj based
/// on this object
OsiBranchingObject *CouenneComplObject::createBranch (OsiSolverInterface *solver, 
						      const OsiBranchingInformation *info, 
						      int way) const {

  expression **args = reference_ -> Image () -> ArgList ();

  /*  printf ("creating CCobj: %d %d.%d\n", reference_ -> Index (), 
	  args [0] -> Index (),
	  args [1] -> Index ());*/

  return new CouenneComplBranchingObject (solver, this, jnlst_,
					  args [0],
					  args [1], 
					  way, 0, doFBBT_, doConvCuts_);
}

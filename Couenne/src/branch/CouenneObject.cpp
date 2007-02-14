/*
 * Name:    CouenneObject.cpp
 * Author:  Pietro Belotti
 * Purpose: Base object for variables (to be used in branching)
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <CouenneObject.hpp>
#include <CouenneBranchingObject.hpp>

/// return something nonnegative (better in [0,1]). 0 means feasible, 
double CouenneObject::infeasibility (const OsiBranchingInformation*, int &) const
  {return fabs ((*reference_) () - (*(reference_ -> Image ())) ());}

/// fix integer coordinates of current integer feasible solution
double CouenneObject::feasibleRegion (OsiSolverInterface *solver, 
				      const OsiBranchingInformation *info) const {

  // fix integer coordinate of current integer feasible solution
  int index = reference_ -> getFixIndex ();
  solver -> setColLower (index, info -> solution_ [index]);
  solver -> setColUpper (index, info -> solution_ [index]);

  return 0.;
}

/// apply the branching rule
OsiBranchingObject* CouenneObject::createBranch (OsiSolverInterface *, 
						 const OsiBranchingInformation *, int) const {

  return new CouenneBranchingObject (reference_);
}

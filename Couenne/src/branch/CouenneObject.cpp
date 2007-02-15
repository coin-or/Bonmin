/*
 * Name:    CouenneObject.cpp
 * Author:  Pietro Belotti
 * Purpose: Base object for variables (to be used in branching)
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <CouenneObject.hpp>
#include <CouenneBranchingObject.hpp>


/// return difference between current value
double CouenneObject::infeasibility (const OsiBranchingInformation*, int &) const
  {return fabs ((*reference_) () - (*(reference_ -> Image ())) ());}


/// fix integer coordinates of current integer feasible solution
double CouenneObject::feasibleRegion (OsiSolverInterface *solver, 
				      const OsiBranchingInformation *info) const {

  // 
  int index = reference_ -> getFixVar () -> Index ();

  // 
  solver -> setColLower (index, info -> solution_ [index]);
  solver -> setColUpper (index, info -> solution_ [index]);

  return 0.;
}


/// apply the branching rule
OsiBranchingObject* CouenneObject::createBranch (OsiSolverInterface *, 
						 const OsiBranchingInformation *, 
						 int) const {

  return new CouenneBranchingObject (reference_ -> getFixVar ());
}

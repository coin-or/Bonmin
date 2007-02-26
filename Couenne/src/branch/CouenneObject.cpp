/*
 * Name:    CouenneObject.cpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Base object for variables (to be used in branching)
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <CouenneObject.hpp>
#include <CouenneBranchingObject.hpp>


/// return difference between current value
double CouenneObject::infeasibility (const OsiBranchingInformation *info, int &) const {

  int index = reference_ -> Image () -> getFixVar () -> Index ();

  // if branched-upon variable has a narrow interval, it is not worthy
  // to branch on it

  if ((fabs (info -> lower_    [index] - info -> upper_    [index]) < COUENNE_EPS)
      //|| (fabs (info -> solution_ [index] - info -> upper_    [index]) < COUENNE_EPS)
      //|| (fabs (info -> solution_ [index] - info -> lower_    [index]) < COUENNE_EPS)
      )
    return 0.;

  const double & expr = (*(reference_ -> Image ())) ();
  const double & var  = expression::Variable (reference_ -> Index ());

  bool verbose = 0;

  if (verbose) {

    reference_ -> print(std::cout);
    std::cout<<" = ";
    reference_->Image()->print(std::cout);
    std::cout<<std::endl;
    std::cout<<expr<<" =(?) "<<var<<std::endl;
  }

  // otherwise, return real value of difference w - f(x)
  return fabs (var - expr);
}


/// fix integer coordinates of current integer feasible solution
double CouenneObject::feasibleRegion (OsiSolverInterface *solver, 
				      const OsiBranchingInformation *info) const {

  // get current value of the branching variable
  int    index = reference_ -> getFixVar () -> Index ();
  double val   = info -> solution_ [index];

  // fix that variable to its current value
  solver -> setColLower (index, val);
  solver -> setColUpper (index, val);

  return 0.;
}


/// apply the branching rule
OsiBranchingObject* CouenneObject::createBranch (OsiSolverInterface *, 
						 const OsiBranchingInformation *, 
						 int) const {

  return new CouenneBranchingObject (reference_ -> Image () -> getFixVar ());
}

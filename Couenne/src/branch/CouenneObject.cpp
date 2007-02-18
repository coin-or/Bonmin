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
double CouenneObject::infeasibility (const OsiBranchingInformation*, int &) const
  {
    const double & expr = (*(reference_ -> Image ())) ();
    const double & var = expression::Variable (reference_ -> Index ());
    bool verbose = 0;
    if(verbose){
      reference_->print(std::cout);
      std::cout<<" = ";
      reference_->Image()->print(std::cout);
      std::cout<<std::endl;
      std::cout<<expr<<" =(?) "<<var<<std::endl;
    }
    return fabs (var - expr);}


/// fix integer coordinates of current integer feasible solution
double CouenneObject::feasibleRegion (OsiSolverInterface *solver, 
				      const OsiBranchingInformation *info) const {

  // 
  int    index = reference_ -> getFixVar () -> Index ();
  double val   = info -> solution_ [index];

  // 
  solver -> setColLower (index, val);
  solver -> setColUpper (index, val);

  return 0.;
}


/// apply the branching rule
OsiBranchingObject* CouenneObject::createBranch (OsiSolverInterface *, 
						 const OsiBranchingInformation *, 
						 int) const {
  /*  printf ("create branch for aux x%d (branch on x%d)\n", 
	  reference_ -> Index(), 
	  reference_ -> Image() -> getFixVar () ->Index());
  */
  return new CouenneBranchingObject (reference_ -> Image () -> getFixVar ());
}

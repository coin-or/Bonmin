/*
 * Name:    CouenneObject.cpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Base object for variables (to be used in branching)
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <CouenneObject.hpp>
#include <exprGroup.h>
#include <CouenneBranchingObject.hpp>


/// return difference between current value
double CouenneObject::infeasibility (const OsiBranchingInformation *info, int &) const {

  int index = reference_ -> Image () -> getFixVar () -> Index ();

  //  if (index < 0) return 0;

  //printf("vars: "); for (int i=0;i<18;i++) printf("%+7.1f ",expression::Variable(i)); printf ("\n");

  expression::update (const_cast <CouNumber *> (info -> solution_),
		      const_cast <CouNumber *> (info -> lower_),
		      const_cast <CouNumber *> (info -> upper_));

  //printf("info: "); for (int i=0;i<18;i++) printf("%+7.1f ",expression::Variable(i)); printf ("\n");

  // if branched-upon variable has a narrow interval, it is not worth
  // to branch on it

  /*
  if ((fabs (info -> lower_ [index] - 
	     info -> upper_ [index]) < COUENNE_EPS)
      //|| (fabs (info -> solution_ [index] - info -> upper_    [index]) < COUENNE_EPS)
      //|| (fabs (info -> solution_ [index] - info -> lower_    [index]) < COUENNE_EPS)
      )
    return 0.;
  */

  const double & expr = (*(reference_ -> Image ())) (), 
               & var  = (*reference_) ();
    //expression::Variable (reference_ -> Index ());

  if (0) {

    reference_             -> print (std::cout); std::cout << " = ";
    reference_ -> Image () -> print (std::cout);

    printf (". Inf: = |%.2f - %.2f| [%.2f,%.2f]= %.2f\n", 
	    var, expr, 
	    expression::Lbound (reference_ -> Index ()),
	    expression::Ubound (reference_ -> Index ()),
	    fabs (var - expr));
  }

  CouNumber delta = fabs (var - expr);

  /// avoid branching on very small deltas
  if (delta < COUENNE_EPS) delta = 0.;

  // otherwise, return real value of difference w - f(x)
  return delta;
}


/// fix integer coordinates of current integer feasible solution
double CouenneObject::feasibleRegion (OsiSolverInterface *solver, 
				      const OsiBranchingInformation *info) const {
  int index = reference_ -> Index ();

  // should never happen...
  if (index < 0) return 0;

  double val = info -> solution_ [index];

  // fix that variable to its current value

  solver -> setColLower (index, val);
  solver -> setColUpper (index, val);

  expression * expr = reference_ -> Image();

  if (expr -> Argument ()){ // unary function

    index = expr -> Argument () -> Index ();

    if (index > -1) {

      val = info -> solution_ [index];

      solver -> setColLower (index, val);
      solver -> setColUpper (index, val);
    }
  }
  else // n-ary function
    if (expr -> ArgList ()) {

      expression ** args = expr -> ArgList ();
      int nargs = expr -> nArgs();

      for (int i = 0 ; i < nargs ; i++) {

	index = args [i] -> Index();

	if (index > -1) {

	  val = info -> solution_ [index];
	  solver -> setColLower (index, val);
	  solver -> setColUpper (index, val);
	}
      }
    }

  // last case: exprGroup, the linear terms are not handled

  if (expr -> code () == COU_EXPRGROUP) {

    exprGroup *e = dynamic_cast <exprGroup *> (expr);
    int *indices = e -> getIndices (), index;

    for (int i=0; (index = indices [i]) >= 0; i++) {
      
      val = info -> solution_ [index];

      solver -> setColLower (index, val);
      solver -> setColUpper (index, val);
    }
  }

  return 0.;
}


/// apply the branching rule
OsiBranchingObject* CouenneObject::createBranch (OsiSolverInterface *si, 
						 const OsiBranchingInformation *, 
						 int) const {
  if (0) {
    printf ("CO::createBranch: ");
    reference_ -> print (std::cout);
    printf (" = ");
    reference_ -> Image () -> print (std::cout);
    printf (" --> branch on ");
    reference_ -> Image () -> getFixVar () -> print (std::cout);
    printf ("\n");
  }

  return new CouenneBranchingObject (reference_ -> Image () -> getFixVar ());
}

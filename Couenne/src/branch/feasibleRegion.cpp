/*
 * Name:    feasibleRegion.cpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Implement feasibleRegion() method of CouenneObject class
 *
 * (C) Carnegie-Mellon University, 2006-08.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"

#include "CouenneSolverInterface.hpp"
#include "CouenneProblem.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"
#include "CouenneThreeWayBranchObj.hpp"

#include "exprGroup.hpp"
#include "exprQuad.hpp"
#include "lqelems.hpp"

#define TOL 0.

/// fix integer coordinates of current integer feasible solution
double CouenneObject::feasibleRegion (OsiSolverInterface *solver, 
				      const OsiBranchingInformation *info) const {
  int index = reference_ -> Index ();

  assert (index >= 0);

  double val = info -> solution_ [index];

  // fix that variable to its current value
  solver -> setColLower (index, val-TOL);
  solver -> setColUpper (index, val+TOL);

  expression *expr = reference_ -> Image ();

  if (!expr) return 0.;

  // fix all variables upon which this auxiliary depends

  // expr is surely nonlinear, so it's useless to check if it is an
  // exprAux, w1:=w0

  if (expr -> Type () == UNARY) { // unary function

    index = expr -> Argument () -> Index ();

    if (index >= 0) {
      val = info -> solution_ [index];
      solver -> setColLower (index, val-TOL);
      solver -> setColUpper (index, val+TOL);
    }
  }
  else // n-ary function

    if (expr -> Type () == N_ARY) {

      expression ** args = expr -> ArgList ();
      int nargs = expr -> nArgs ();

      for (register int i=0; i < nargs; i++) {

	if ((index = args [i] -> Index()) >= 0) {
	  val = info -> solution_ [index];
	  solver -> setColLower (index, val-TOL);
	  solver -> setColUpper (index, val+TOL);
	}
      }
    }

  // last cases: exprGroup/Quad, must handle the linear/quadratic terms
  if ((expr -> code () == COU_EXPRGROUP) ||
      (expr -> code () == COU_EXPRQUAD)) {

    exprGroup *e = dynamic_cast <exprGroup *> (expr);

    exprGroup::lincoeff &lcoe = e -> lcoeff ();

    for (exprGroup::lincoeff::iterator el = lcoe.begin (); el != lcoe.end (); ++el) {
      int index = el -> first -> Index ();
      val = info -> solution_ [index];
      solver -> setColLower (index, val-TOL);
      solver -> setColUpper (index, val+TOL);
    }

    // take care of quadratic terms
    if (expr -> code () == COU_EXPRQUAD) {

      exprQuad *e = dynamic_cast <exprQuad *> (expr);

      exprQuad::sparseQ q = e -> getQ ();

      for (exprQuad::sparseQ::iterator row = q.begin (); 
	   row != q.end (); ++row) {

	int xind = row -> first -> Index ();

	val = info -> solution_ [xind];
	solver -> setColLower (xind, val-TOL);
	solver -> setColUpper (xind, val+TOL);

	for (exprQuad::sparseQcol::iterator col = row -> second.begin ();
	     col != row -> second.end (); ++col) {

	  int yind = col -> first -> Index ();

	  val = info -> solution_ [yind];
	  solver -> setColLower (yind, val-TOL);
	  solver -> setColUpper (yind, val+TOL);
	}
      }
    }
  }

  return 0.;
}

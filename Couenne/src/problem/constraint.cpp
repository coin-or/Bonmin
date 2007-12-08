/*
 * Name:    constraint.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the classes CouenneConstraint and LinearConstraint
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneTypes.hpp"
#include "CouenneProblemElem.hpp"


// output nonlinear constraint

void CouenneConstraint::print (std::ostream &out) {

  bool samebounds = ((lb_ -> Type () == CONST) &&
		     (ub_ -> Type () == CONST) && 
		     (fabs (lb_ -> Value () - ub_ -> Value ()) < COUENNE_EPS));

  // left hand side (a in a <= h(x) <= b)

  if (lb_ && 
      !samebounds &&
      ((lb_ -> Type  () != CONST) ||
       (lb_ -> Value () > - COUENNE_INFINITY))) {

    lb_ -> print (out); fflush (stdout);
    out  << " <= "; fflush (stdout);
  }

  // body: h(x) in a <= h(x) <= b

  body_ -> print (out); fflush (stdout);

  // right hand side

  if (ub_ && ((ub_ -> Type  () != CONST) || 
	      (ub_ -> Value () <  COUENNE_INFINITY))) {

    out << ' ';
    if (!samebounds) out << "<";
    out << "= "; fflush (stdout);
    ub_ -> print (out); fflush (stdout);
  } 

  out << std::endl;
}

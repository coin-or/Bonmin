/*
 * Name:    constraint.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the classes CouenneConstraint and LinearConstraint
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>
#include <map>

#include <CouenneTypes.h>
#include <CouenneProblem.hpp>
#include <CouenneProblemElem.hpp>


// output nonlinear constraint

void CouenneConstraint::print (std::ostream &out = std::cout) {

  bool samebounds = ((lb_ -> Type () == CONST) &&
		     (ub_ -> Type () == CONST) && 
		     (fabs (lb_ -> Value () - ub_ -> Value ()) < COUENNE_EPS));

  if (lb_ && 
      !samebounds &&
      ((lb_ -> Type  () != CONST) ||
       (lb_ -> Value () > - COUENNE_INFINITY))) {

    lb_ -> print (out); fflush (stdout);
    out  << " <= "; fflush (stdout);
  }

  body_ -> print (out); fflush (stdout);

  if (ub_ && ((ub_ -> Type  () != CONST) || 
	      (ub_ -> Value () <  COUENNE_INFINITY))) {

    out << ' ';
    if (!samebounds) out << "<";
    out << "= "; fflush (stdout);
    ub_ -> print (out); fflush (stdout);
  } 

  out << std::endl;
}

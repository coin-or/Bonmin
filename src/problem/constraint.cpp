/*
 * Name:    constraint.C
 * Author:  Pietro Belotti
 * Purpose: methods of the classes CouenneConstraint and LinearConstraint
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>
#include <map>

#include <CouenneTypes.h>
#include <CouenneProblem.h>
#include <CouenneProblemElem.h>


// output nonlinear constraint

void CouenneConstraint::print (std::ostream &out = std::cout) {

  if (lb_ && ((lb_ -> Type  () != CONST) ||
	      (lb_ -> Value () > - COUENNE_INFINITY))) {

    lb_ -> print (out); fflush (stdout);
    out  << " <= "; fflush (stdout);
  }

  body_ -> print (out); fflush (stdout);

  if (ub_ && ((ub_ -> Type  () != CONST) || 
	      (ub_ -> Value () <  COUENNE_INFINITY))) {

    out << " <= "; fflush (stdout);
    ub_ -> print (out); fflush (stdout);
  } 

  out << std::endl;
}


// output linear constraint

void LinearConstraint::print (std::ostream &out = std::cout) {
  /*
  if (ub_ && lb_ &&
      ((lb_ -> Type  () != CONST) ||
       (lb_ -> Value () > - COUENNE_INFINITY))) {

    lb_ -> print (out); 
    out << " <= ";
  }
  */

  for (int i = 0; i<nterms_; i++) {

    if (i) out << " +";
    coeff_ [i] -> print (out); fflush (stdout);
    out << "*x_" << indices_ [i];
  }

  if      (sign_ == COUENNE_EQ) out << " = "; 
  else if (sign_ == COUENNE_GE) out << " >= ";
  else                          out << " <= ";

  rhs_ -> print (out);

  /*
  if (ub_ && ((1) || (ub_ -> Type  () != CONST) || 
	      (ub_ -> Value () <  COUENNE_INFINITY)))
    ub_ -> print (out);
  else
    if  (lb_ && ((1) || (lb_ -> Type  () != CONST) ||
		 (lb_ -> Value () > - COUENNE_INFINITY)))
      lb_ -> print (out); 
    else {
      printf (" << ");
      if (lb_) lb_ -> print (out); 
      else printf ("(null!)");
      printf (" -- ");
      if (ub_) ub_ -> print (out); 
      else printf ("(null!)");
      printf (" >> ");
    }
  */

  out << std::endl;
}

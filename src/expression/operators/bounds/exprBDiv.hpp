/*
 * Name:    exprBDiv.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of operators to compute lower/upper bounds of divisions
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRBDIV_H
#define COUENNE_EXPRBDIV_H

#include "exprOp.hpp"


/// division that avoids NaN's and considers a sign when returning infinity
static inline CouNumber safeDiv (register CouNumber a, register CouNumber b, int sign) {

  if (fabs (a) < COUENNE_EPS) return 0;
  //    if (fabs (b) < COUENNE_EPS)) return 0;
  //    else return 0 

  if (fabs (b) < COUENNE_EPS) return ((sign < 0) ? -COUENNE_INFINITY :  COUENNE_INFINITY);

  if (a >  COUENNE_INFINITY) return ((sign < 0) ? -COUENNE_INFINITY :  COUENNE_INFINITY);
  if (a < -COUENNE_INFINITY) return ((sign < 0) ? -COUENNE_INFINITY :  COUENNE_INFINITY); 

  return a/b;
}


///  class to compute lower bound of a fraction based on the bounds of
///  both numerator and denominator

class exprLBDiv: public exprOp {

 public:

  /// Constructors, destructor
  exprLBDiv  (expression **al, int n): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  /// cloning method
  expression *clone (Domain *d = NULL) const
    {return new exprLBDiv (clonearglist (d), nargs_);}

  /// function for the evaluation of the expression
  CouNumber operator () ();

  /// print position (PRE, INSIDE, POST)
  enum pos printPos () const
    {return PRE;}

  /// print operator
  std::string printOp () const
    {return "LB_Div";}
};


/// compute sum

inline CouNumber exprLBDiv::operator () () {

  register CouNumber n = (*(arglist_ [0])) ();
  register CouNumber N = (*(arglist_ [1])) ();
  register CouNumber d = (*(arglist_ [2])) ();
  register CouNumber D = (*(arglist_ [3])) ();
                                               // (n,N,d,D)     lb 
  if (d > 0)                                   // (?,?,+,+)
    if   (n > 0)    return safeDiv (n,D,-1);      // (+,+,+,+) --> n/D
    else            return safeDiv (n,d,-1);      // (-,?,+,+) --> n/d
  else { // d <= 0
    if      (D > 0) return - COUENNE_INFINITY; // (?,?,-,+) --> unbounded
    else if (N > 0) return safeDiv (N,D,-1);      // (?,+,-,-) --> N/D
    else            return safeDiv (N,d,-1);      // (-,-,-,-) --> N/d
  }
}


/// class to compute upper bound of a fraction based on the bounds of
/// both numerator and denominator

class exprUBDiv: public exprOp {

 public:

  /// Constructors, destructor
  exprUBDiv  (expression **al, int n): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  /// cloning method
  expression *clone (Domain *d = NULL) const
  {return new exprUBDiv (clonearglist (d), nargs_);}

  /// function for the evaluation of the expression
  CouNumber operator () ();

  /// print position (PRE, INSIDE, POST)
  enum pos printPos () const
    {return PRE;}

  /// print operator
  std::string printOp () const
    {return "UB_Div";}
};


/// compute sum

inline CouNumber exprUBDiv::operator () () {

  register CouNumber n = (*(arglist_ [0])) ();
  register CouNumber N = (*(arglist_ [1])) ();
  register CouNumber d = (*(arglist_ [2])) ();
  register CouNumber D = (*(arglist_ [3])) ();

  if (d > 0)                                     // (n,N,d,D)     lb 
    if   (N < 0) return safeDiv (N,D,1);         // (-,-,+,+) --> N/D
    else         return safeDiv (N,d,1);         // (?,+,+,+) --> N/d
  else { // d <= 0
    if      (D > 0) return + COUENNE_INFINITY;   // (?,?,-,+) --> unbounded
    else if (n < 0) return safeDiv (n,D,1);      // (-,?,-,-) --> n/D
    else            return safeDiv (n,d,1);      // (+,+,-,-) --> n/d
  }
}

#endif

/*
 * Name:    exprBDiv.h
 * Author:  Pietro Belotti
 * Purpose: definition of operators to compute lower/upper bounds of divisions
 *
 * (C) Pietro Belotti 2006. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRBDIV_H
#define COUENNE_EXPRBDIV_H

#include <exprOp.h>


/// division that avoids NaN's 
inline CouNumber safeDiv (register CouNumber a, register CouNumber b) {

  if ((fabs (a) < 1e-10) && (fabs (b) < 1e-10)) return 0;
  if (fabs (a) > 1e20) {
    if (a > 1e20) return ((b < 0) ? -COUENNE_INFINITY :  COUENNE_INFINITY);
    else          return ((b < 0) ?  COUENNE_INFINITY : -COUENNE_INFINITY);
  }
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
  expression *clone () const
    {return new exprLBDiv (clonearglist (), nargs_);}

  /// function for the evaluation of the expression
  CouNumber operator () ();

  /// I/O
  void print (std::ostream &) const;
};


/// output
inline void exprLBDiv::print (std::ostream &out = std::cout) const
{exprOp::print (out, "LB_div", PRE);}


/// compute sum

inline CouNumber exprLBDiv::operator () () {

  //  exprOp:: operator () ();

  register CouNumber n = (*(arglist_ [0])) ();
  register CouNumber N = (*(arglist_ [1])) ();
  register CouNumber d = (*(arglist_ [2])) ();
  register CouNumber D = (*(arglist_ [3])) ();
  /*
  register CouNumber D = *sp--;
  register CouNumber d = *sp--;
  register CouNumber N = *sp--;
  register CouNumber n = *sp--;
  */

  //  printf ("lbdiv: %e %e %e %e\n", n,N,d,D);
                                                     // (n,N,d,D)     lb 
  if (d > COUENNE_EPS) {                             // (?,?,+,+)
    if   (n > 0) return safeDiv (n,D);               // (+,+,+,+) --> n/D
    else         return safeDiv (n,d);               // (-,?,+,+) --> n/d
  } else  // d <= 0
    if      (D > COUENNE_EPS) return - COUENNE_INFINITY; // (?,?,-,+) --> unbounded
    else if (N > COUENNE_EPS) return safeDiv (N,D);      // (?,+,-,-) --> N/D
    else                      return safeDiv (N,d);      // (-,-,-,-) --> N/d
}


///  class to compute lower bound of a fraction based on the bounds of
///  both numerator and denominator

class exprUBDiv: public exprOp {

 public:

  /// Constructors, destructor
  exprUBDiv  (expression **al, int n): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  /// cloning method
  expression *clone () const
    {return new exprUBDiv (clonearglist (), nargs_);}

  /// function for the evaluation of the expression
  CouNumber operator () ();

  /// output
  void print (std::ostream &) const;
};


/// compute sum

inline CouNumber exprUBDiv::operator () () {

  //  exprOp:: operator () ();

  register CouNumber n = (*(arglist_ [0])) ();
  register CouNumber N = (*(arglist_ [1])) ();
  register CouNumber d = (*(arglist_ [2])) ();
  register CouNumber D = (*(arglist_ [3])) ();
  /*
  register CouNumber D = *sp--;
  register CouNumber d = *sp--;
  register CouNumber N = *sp--;
  register CouNumber n = *sp--;
  */

  //  printf ("ubdiv: %e %e %e %e\n", n,N,d,D);
                                                       // (n,N,d,D)     lb 
  if (d > COUENNE_EPS) {                                                     
    if   (N < 0) return safeDiv (N,D);                 // (-,-,+,+) --> N/D
    else         return safeDiv (N,d);                 // (?,+,+,+) --> N/d
  } else { // d <= 0
    if      (D >   COUENNE_EPS) return + COUENNE_INFINITY; // (?,?,-,+) --> unbounded
    else if (n < - COUENNE_EPS) return safeDiv (n,D);  // (-,?,-,-) --> n/D
    else                        return safeDiv (n,d);  // (+,+,-,-) --> n/d
  }
}

/// output
void exprUBDiv::print (std::ostream &out = std::cout) const
{exprOp::print (out, "UB_div", PRE);}

#endif

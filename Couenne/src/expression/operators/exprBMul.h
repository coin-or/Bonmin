/*
 * Name:    exprBMul.h
 * Author:  Pietro Belotti
 * Purpose: definition of operators to compute lower/upper bounds of multiplications
 *
 * (C) Pietro Belotti 2006. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRBMUL_H
#define COUENNE_EXPRBMUL_H

#include <exprOp.h>


// product that avoids NaN's 
inline CouNumber safeProd (register CouNumber a, register CouNumber b) {

  if ((fabs (a) < 1e-10) || (fabs (b) < 1e-10)) 
    return 0.;
  return a*b;
}


//  class to compute lower bound of a fraction based on the bounds of
//  both numerator and denominator

class exprLBMul: public exprOp {

 public:

  // Constructors, destructor
  exprLBMul  (expression **al, int n): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  // cloning method
  expression *clone () const
    {return new exprLBMul (clonearglist (), nargs_);}

  // function for the evaluation of the expression
  CouNumber operator () ();

  // output
  void print (std::ostream &) const;
};


// compute sum

inline CouNumber exprLBMul::operator () () {

  register CouNumber n = (*(arglist_ [0])) ();
  register CouNumber N = (*(arglist_ [1])) ();
  register CouNumber d = (*(arglist_ [2])) ();
  register CouNumber D = (*(arglist_ [3])) ();

  if (d>=0)
    if   (n>=0) return safeProd (n,d);
    else        return safeProd (n,D);
  else // d <= 0
    if (N>0) {
      CouNumber Nd = safeProd (N,d), nD;
      if (n<0 && D>0 && 
	  (Nd > (nD = safeProd (n,D)))) return nD;
      else                              return Nd;
    }
    else 
      if (D>0) return safeProd (n,D);
      else     return safeProd (N,D);
}


// output
inline void exprLBMul::print (std::ostream &out = std::cout) const
{exprOp::print (out, "LB_mul", PRE);}


//  class to compute lower bound of a fraction based on the bounds of
//  both numerator and denominator

class exprUBMul: public exprOp {

 public:

  // Constructors, destructor
  exprUBMul  (expression **al, int n): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  // cloning method
  expression *clone () const
    {return new exprUBMul (clonearglist (), nargs_);}

  // function for the evaluation of the expression
  CouNumber operator () ();

  // output
  void print (std::ostream &) const;
};



// compute sum

inline CouNumber exprUBMul::operator () () {

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

  if (d>0)
    if (N<0) return safeProd (N,d);
    else     return safeProd (N,D);
  else // d <= 0
    if (n<0) {
      CouNumber nd = safeProd (n,d), ND;
      if (N>0 && D>0 && 
	  ((ND = safeProd (N,D)) > nd)) return ND;
      else                              return nd;
    }
    else 
      if (D>0) return safeProd (N,D);
      else     return safeProd (n,D);
}


// output

void exprUBMul::print (std::ostream &out = std::cout) const
{exprOp::print (out, "UB_mul", PRE);}

#endif

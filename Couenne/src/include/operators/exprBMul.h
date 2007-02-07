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


//  class to compute lower bound of a fraction based on the bounds of
//  both numerator and denominator

class exprLBMul: public exprOp {

 public:

  // Constructors, destructor
  exprLBMul  (expression **al, int n): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  ~exprLBMul () {}

  // cloning method
  expression *clone () const
    {return new exprLBMul (clonearglist (), nargs_);}

  // function for the evaluation of the expression
  CouNumber operator () ();

  // String equivalent (for comparisons)
  const std::string name() const {return "LB_mul" + exprOp::name();}

  // output
  void print (std::ostream &) const;
};


// compute sum

inline CouNumber exprLBMul::operator () () {

  exprOp:: operator () ();

  register CouNumber D = *sp--;
  register CouNumber d = *sp--;
  register CouNumber N = *sp--;
  register CouNumber n = *sp--;
  CouNumber nD;

  if (d>0)
    if   (n>0) return n*d;
    else       return n*D;
  else // d <= 0
    if (N>0)
      if (n<0 && D>0 && N*d > (nD = n*D)) return nD;
      else                                return N*d;
    else 
      if (D>0) return n*D;
      else     return N*D;
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

  ~exprUBMul () {}

  // cloning method
  expression *clone () const
    {return new exprUBMul (clonearglist (), nargs_);}

  // function for the evaluation of the expression
  CouNumber operator () ();

  // String equivalent (for comparisons)
  std::string name() {return "UB_mul" + exprOp::name();}

  // output
  void print (std::ostream &) const;
};



// compute sum

inline CouNumber exprUBMul::operator () () {

  exprOp:: operator () ();

  register CouNumber D = *sp--;
  register CouNumber d = *sp--;
  register CouNumber N = *sp--;
  register CouNumber n = *sp--;
  CouNumber ND;

  if (d>0)
    if (N<0) return N*d;
    else     return N*D;
  else // d <= 0
    if (n<0) 
      if (N>0 && D>0 && ((ND=N*D) > n*d)) return ND;
      else                                return n*d;
    else 
      if (D>0) return N*D;
      else     return n*D;
}


// output

inline void exprUBMul::print (std::ostream &out = std::cout) const
{exprOp::print (out, "UB_mul", PRE);}

#endif

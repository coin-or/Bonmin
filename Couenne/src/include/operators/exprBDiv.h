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


//  class to compute lower bound of a fraction based on the bounds of
//  both numerator and denominator

class exprLBDiv: public exprOp {

 public:

  // Constructors, destructor
  exprLBDiv  (expression **al, int n): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  ~exprLBDiv () {}

  // cloning method
  expression *clone () const
    {return new exprLBDiv (clonearglist (), nargs_);}

  // function for the evaluation of the expression
  CouNumber operator () ();

  // output
  void print (std::ostream &);
};


// output
inline void exprLBDiv::print (std::ostream &out = std::cout)
{exprOp::print (out, "LB_div", PRE);}


// compute sum

inline CouNumber exprLBDiv::operator () () {

  exprOp:: operator () ();

  register CouNumber D = *sp--;
  register CouNumber d = *sp--;
  register CouNumber N = *sp--;
  register CouNumber n = *sp--;
                                                     // (n,N,d,D)     lb 
  if (d > COUENNE_EPS) {                             // (?,?,+,+)
    if   (n > 0) return n/D;                         // (+,+,+,+) --> n/D
    else         return n/d;                         // (-,?,+,+) --> n/d
  } else  // d <= 0
    if      (D > COUENNE_EPS) return - COUENNE_INFINITY; // (?,?,-,+) --> unbounded
    else if (N > COUENNE_EPS) return N/D;                // (?,+,-,-) --> N/D
    else                      return N/d;                // (-,-,-,-) --> N/d
}


//  class to compute lower bound of a fraction based on the bounds of
//  both numerator and denominator

class exprUBDiv: public exprOp {

 public:

  // Constructors, destructor
  exprUBDiv  (expression **al, int n): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  ~exprUBDiv () {}

  // cloning method
  expression *clone () const
    {return new exprUBDiv (clonearglist (), nargs_);}

  // function for the evaluation of the expression
  CouNumber operator () ();

  // output
  void print (std::ostream &);
};


// compute sum

inline CouNumber exprUBDiv::operator () () {

  exprOp:: operator () ();

  register CouNumber D = *sp--;
  register CouNumber d = *sp--;
  register CouNumber N = *sp--;
  register CouNumber n = *sp--;
                                                       // (n,N,d,D)     lb 
  if (d > COUENNE_EPS) {                                                     
    if   (N < 0) return N/D;                           // (-,-,+,+) --> N/D
    else         return N/d;                           // (?,+,+,+) --> N/d
  } else { // d <= 0
    if      (D >   COUENNE_EPS) return + COUENNE_INFINITY; // (?,?,-,+) --> unbounded
    else if (n < - COUENNE_EPS) return n/D;              // (-,?,-,-) --> n/D
    else                        return n/d;              // (+,+,-,-) --> n/d
  }
}

// output
inline void exprUBDiv::print (std::ostream &out = std::cout)
{exprOp::print (out, "UB_div", PRE);}

#endif

/*
 * Name:    exprBCos.h
 * Author:  Pietro Belotti
 * Purpose: definition of operators to compute lower/upper bounds of cosines
 *
 * (C) Pietro Belotti 2006. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRBCOS_H
#define COUENNE_EXPRBCOS_H

#include <exprOp.h>
#include <math.h>


//  class to compute lower bound of a cosine based on the bounds of
//  its arguments

class exprLBCos: public exprOp {

 public:

  // Constructors, destructor
  exprLBCos (expression *lb, expression *ub): 
    exprOp (new expression * [2], 2) {
    arglist_ [0] = lb;
    arglist_ [1] = ub;
  } //< non-leaf expression, with argument list

  // cloning method
  expression *clone () const
    {return new exprLBCos (arglist_ [0] -> clone (), 
			   arglist_ [1] -> clone ());}

  // function for the evaluation of the expression
  CouNumber operator () ();

  // output
  void print (std::ostream &out = std::cout) 
    {exprOp::print (out, "LB_Cos", PRE);}
};


// compute sum

inline CouNumber exprLBCos::operator () () {

  register CouNumber l = (*(arglist_ [0])) ();
  register CouNumber u = (*(arglist_ [1])) ();

  if ((u - l > 2 * M_PI) ||      // 1) interval spans whole cycle
      (floor (l/2/M_PI - 0.5) < // 2) there is a 3/2 pi + 2k pi in between
       floor (u/2/M_PI - 0.5))) 
    return -1;

  return mymin (sin (l), sin (u));
}

///////////////////////////////////////////////////////////////////////////////

//  class to compute lower bound of a cosine based on the bounds of
//  its arguments

class exprUBCos: public exprOp {

 public:

  // Constructors, destructor
  exprUBCos (expression *lb, expression *ub): 
    exprOp (new expression * [2], 2) {
    arglist_ [0] = lb;
    arglist_ [1] = ub;
  } //< non-leaf expression, with argument list

  // cloning method
  expression *clone () const
    {return new exprUBCos (arglist_ [0] -> clone (), 
			   arglist_ [1] -> clone ());}

  // function for the evaluation of the expression
  CouNumber operator () ();

  // output
  void print (std::ostream &out = std::cout) 
    {exprOp::print (out, "UB_Cos", PRE);}
};


// compute sum

inline CouNumber exprUBCos::operator () () {

  register CouNumber l = (*(arglist_ [0])) ();
  register CouNumber u = (*(arglist_ [1])) ();

  if ((u - l > 2 * M_PI) || // 1) interval spans whole cycle
      (floor (l/2/M_PI) <   // 2) there is a 3/2 pi + 2k pi in between
       floor (u/2/M_PI))) 
    return 1;

  return mymax (sin (l), sin (u));
}

#endif

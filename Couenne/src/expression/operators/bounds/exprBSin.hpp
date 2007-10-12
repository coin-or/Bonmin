/*
 * Name:    exprBSin.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of operators to compute lower/upper bounds of sines
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRBSIN_HPP
#define COUENNE_EXPRBSIN_HPP

#include <exprOp.hpp>

#include <math.h>
#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

//  class to compute lower bound of a sine based on the bounds on its
//  arguments

class exprLBSin: public exprOp {

 public:

  // Constructors, destructor
  exprLBSin (expression *lb, expression *ub): 
    exprOp (new expression * [2], 2) {
    arglist_ [0] = lb;
    arglist_ [1] = ub;
  } //< non-leaf expression, with argument list

  // cloning method
  expression *clone () const
    {return new exprLBSin  (arglist_ [0] -> clone (), 
			    arglist_ [1] -> clone ());}

  // function for the evaluation of the expression
  CouNumber operator () ();

  /// print operator
  std::string printOp () const
    {return "LB_Sin";}

  /// print position
  enum pos printPos () const
    {return PRE;}
};


// compute sum

inline CouNumber exprLBSin::operator () () {

  register CouNumber 
    l = (*(arglist_ [0])) (),
    u = (*(arglist_ [1])) ();

  CouNumber pi2 = 2 * M_PI;
 
  if ((u - l > pi2) ||        // 1) interval spans whole cycle
      (floor (l/pi2 - 0.75) < // 2) there is a 3/2 pi + 2k pi in between
       floor (u/pi2 - 0.75))) 
    return -1.;

  return mymin (sin (l), sin (u));
}


///////////////////////////////////////////////////////////////////////////////

//  class to compute lower bound of a sine based on the bounds on its
//  arguments

class exprUBSin: public exprOp {

 public:

  // Constructors, destructor
  exprUBSin (expression *lb, expression *ub): 
    exprOp (new expression * [2], 2) {
    arglist_ [0] = lb;
    arglist_ [1] = ub;
  } //< non-leaf expression, with argument list

  // cloning method
  expression *clone () const
    {return new exprUBSin  (arglist_ [0] -> clone (), 
			    arglist_ [1] -> clone ());}

  // function for the evaluation of the expression
  CouNumber operator () ();

  /// print operator
  std::string printOp () const
    {return "UB_Sin";}

  /// print position
  enum pos printPos () const
    {return PRE;}
};


// compute sum

inline CouNumber exprUBSin::operator () () {

  register CouNumber 
    l = (*(arglist_ [0])) (),
    u = (*(arglist_ [1])) ();

  CouNumber pi2 = 2 * M_PI;

  if ((u - l > pi2) ||        // 1) interval spans whole cycle
      (floor (l/pi2 - 0.25) < // 2) there is a pi/2 + 2k pi in between
       floor (u/pi2 - 0.25))) 
    return 1.;

  return mymax (sin (l), sin (u));
}

#endif

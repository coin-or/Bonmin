/*
 * Name:    exprBCos.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of operators to compute lower/upper bounds of cosines
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRBCOS_HPP
#define COUENNE_EXPRBCOS_HPP

#include "exprOp.hpp"
#include "CoinHelperFunctions.hpp"

#include <math.h>
#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif


///  class to compute lower bound of a cosine based on the bounds of
///  its arguments

class exprLBCos: public exprOp {

 public:

  /// Constructors, destructor
  exprLBCos (expression *lb, expression *ub): 
    exprOp (new expression * [2], 2) {
    arglist_ [0] = lb;
    arglist_ [1] = ub;
  } //< non-leaf expression, with argument list

  /// cloning method
  expression *clone (Domain *d = NULL) const
  {return new exprLBCos (arglist_ [0] -> clone (d), 
			 arglist_ [1] -> clone (d));}

  /// function for the evaluation of the expression
  CouNumber operator () ();

  /// print position (PRE, INSIDE, POST)
  enum pos printPos () const
    {return PRE;}

  /// print operator
  std::string printOp () const
    {return "LB_Cos";}
};


/// compute sum

inline CouNumber exprLBCos::operator () () {

  register CouNumber 
    l = (*(arglist_ [0])) (),
    u = (*(arglist_ [1])) ();

  CouNumber pi2 = 2 * M_PI;
 
  if ((u - l > pi2) ||       // 1) interval spans whole cycle
      (floor (l/pi2 - 0.5) < // 2) there is a pi + 2k pi between l and u
       floor (u/pi2 - 0.5))) 
    return -1.;

  return CoinMin (cos (l), cos (u));
}


///////////////////////////////////////////////////////////////////////////////

///  class to compute lower bound of a cosine based on the bounds of
///  its arguments

class exprUBCos: public exprOp {

 public:

  /// Constructors, destructor
  exprUBCos (expression *lb, expression *ub): 
    exprOp (new expression * [2], 2) {
    arglist_ [0] = lb;
    arglist_ [1] = ub;
  } //< non-leaf expression, with argument list

  /// cloning method
  expression *clone (Domain *d = NULL) const
  {return new exprUBCos (arglist_ [0] -> clone (d), 
			 arglist_ [1] -> clone (d));}

  /// function for the evaluation of the expression
  CouNumber operator () ();

  /// print operator
  std::string printOp () const
    {return "UB_Cos";}

  /// print position (PRE, INSIDE, POST)
  enum pos printPos () const
    {return PRE;}
};


/// compute sum

inline CouNumber exprUBCos::operator () () {

  register CouNumber 
    l = (*(arglist_ [0])) (),
    u = (*(arglist_ [1])) ();

  CouNumber pi2 = 2 * M_PI;

  if ((u - l > pi2) || // 1) interval spans whole cycle
      (floor (l/pi2) < // 2) there is a 3/2 pi + 2k pi between l and u
       floor (u/pi2))) 
    return 1.;

  return CoinMax (cos (l), cos (u));
}

#endif

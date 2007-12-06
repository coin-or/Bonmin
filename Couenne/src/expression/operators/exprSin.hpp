/*
 * Name:    exprSin.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the sine of a function
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRSIN_HPP
#define COUENNE_EXPRSIN_HPP

#include <math.h>

#include "exprUnary.hpp"
#include "exprConst.hpp"


/// specify which trigonometric function is dealt with in trigEnvelope
enum cou_trig {COU_SINE, COU_COSINE};


/// normalize angle within [0,b] (typically, pi or 2pi)
inline CouNumber modulo (register CouNumber a, register CouNumber b)
  {return a - b * floor (a/b);}


/// generalized procedure for both sine and cosine
CouNumber trigSelBranch (const CouenneObject *obj, 
			 const OsiBranchingInformation *info,
			 int &ind, 
			 double * &brpts, 
			 int &way,
			 enum cou_trig type);


/// generalized implied bound procedure for sine/cosine
bool trigImpliedBound (enum cou_trig, int, int, CouNumber *, CouNumber *, t_chg_bounds *);


/// class for sin f(x)

class exprSin: public exprUnary {

 public:

  /// Constructors, destructor
  exprSin (expression *al): 
    exprUnary (al) {} //< non-leaf expression, with argument list

  /// cloning method
  expression *clone () const
  {return new exprSin (argument_ -> clone ());}

  //// the operator's function
  inline unary_function F () 
  {return sin;}

  /// print operator
  std::string printOp () const
  {return "sin";}

  /// differentiation
  expression *differentiate (int index); 

  /// Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  /// generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg, 
		     t_chg_bounds * = NULL, int = -1, 
		     CouNumber = -COUENNE_INFINITY, 
		     CouNumber =  COUENNE_INFINITY);

  /// code for comparisons
  virtual enum expr_type code () 
  {return COU_EXPRSIN;}

  /// implied bound processing
  bool impliedBound (int index, CouNumber *l, CouNumber *u, t_chg_bounds *chg)
  {return trigImpliedBound (COU_SINE, index, argument_ -> Index (), l, u, chg);}

  /// Set up branching object by evaluating many branching points for
  /// each expression's arguments
  CouNumber selectBranch (const CouenneObject *obj, const OsiBranchingInformation *info,
			  int &ind, double * &brpts, int &way)
  {return trigSelBranch (obj, info, ind, brpts, way, COU_SINE);}
};

#endif

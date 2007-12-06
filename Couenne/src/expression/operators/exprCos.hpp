/*
 * Name:    exprCos.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of cosine 
 *
 * (C) Carnegie-Mellon University, 2006-07. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRCOS_HPP
#define COUENNE_EXPRCOS_HPP

#include "exprSin.hpp"

/// class cosine

class exprCos: public exprUnary {

 public:

  /// constructor, destructor
  exprCos (expression *al):
    exprUnary (al) {}

  /// cloning method
  expression *clone () const
  {return new exprCos (argument_ -> clone ());}

  //// the operator's function
  inline unary_function F () 
  {return cos;}

  /// print operator
  std::string printOp () const
    {return "cos";}

  /// obtain derivative of expression
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
  {return COU_EXPRCOS;}

  /// Set up branching object by evaluating many branching points for
  /// each expression's arguments
  CouNumber selectBranch (const CouenneObject *obj, const OsiBranchingInformation *info,
			  int &ind, double * &brpts, int &way)
  {return trigSelBranch (obj, info, ind, brpts, way, COU_COSINE);}
};


/// common convexification method used by both cos and sin
CouNumber trigNewton (CouNumber, CouNumber, CouNumber);

#endif

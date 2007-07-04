/*
 * Name:    exprCos.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of cosine 
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRCOS_H
#define COUENNE_EXPRCOS_H

#include <exprUnary.hpp>
#include <exprConst.hpp>
#include <CouennePrecisions.h>

#include <math.h>


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
  inline unary_function F () {return cos;}

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
		     t_chg_bounds * = NULL);

  /// code for comparisons
  virtual enum expr_type code () {return COU_EXPRCOS;}
};


/// common convexification method used by both cos and sin

CouNumber trigNewton (CouNumber, CouNumber, CouNumber);

/// convex envelope for sine/cosine 

#endif

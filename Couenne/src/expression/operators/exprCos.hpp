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
  expression *clone (Domain *d = NULL) const
  {return new exprCos (argument_ -> clone (d));}

  //// the operator's function
  inline unary_function F () 
  {return cos;}

  /// print operator
  std::string printOp () const
  {return "cos";}

  /// return l-2 norm of gradient at given point
  inline CouNumber gradientNorm (const double *x) {
    return (argument_ -> Index () < 0) ? 
      0. : fabs (sin (x [argument_ -> Index ()]));
  }

  /// obtain derivative of expression
  expression *differentiate (int index); 

  /// Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  /// generate equality between *this and *w
  void generateCuts (expression *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg, 
		     t_chg_bounds * = NULL, int = -1, 
		     CouNumber = -COUENNE_INFINITY, 
		     CouNumber =  COUENNE_INFINITY);

  /// code for comparisons
  virtual enum expr_type code ()
  {return COU_EXPRCOS;}

  /// implied bound processing
  bool impliedBound (int index, CouNumber *l, CouNumber *u, t_chg_bounds *chg) {

    bool impl = trigImpliedBound (COU_COSINE, index, argument_ -> Index (), l, u, chg);

    if (impl && argument_ -> isInteger ()) {

      int ind = argument_ -> Index ();
      assert (ind >= 0);
      l [ind] = ceil  (l [ind] - COUENNE_EPS);
      u [ind] = floor (u [ind] + COUENNE_EPS);
    }

    return impl;
  }

  /// Set up branching object by evaluating many branching points for
  /// each expression's arguments
  virtual CouNumber selectBranch (const CouenneObject *obj, 
				  const OsiBranchingInformation *info,
				  expression * &var, 
				  double * &brpts, 
 				  double * &brDist, // distance of current LP
					  	    // point to new convexifications
				  int &way)
  {return trigSelBranch (obj, info, var, brpts, brDist, way, COU_COSINE);}

  /// closest feasible points in function in both directions
  virtual void closestFeasible (expression *varind, expression *vardep,
				CouNumber& left, CouNumber& right) const;

  /// can this expression be further linearized or are we on its
  /// concave ("bad") side
  virtual bool isCuttable (CouenneProblem *problem, int index) const
  {return false;}
};


/// common convexification method used by both cos and sin
CouNumber trigNewton (CouNumber, CouNumber, CouNumber);

#endif

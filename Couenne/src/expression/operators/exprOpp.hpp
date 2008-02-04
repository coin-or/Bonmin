/*
 * Name:    exprOpp.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the opposite -f(x) of a function
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPROPP_HPP
#define COUENNE_EXPROPP_HPP

#include "CouennePrecisions.hpp"
#include "exprUnary.hpp"


/// operator opp: returns the opposite of a number

inline CouNumber opp (register CouNumber arg) 
{return - arg;}


/// class opposite 

class exprOpp: public exprUnary {

 public:

  /// Constructors, destructor
  exprOpp (expression *al): 
    exprUnary (al) {} //< non-leaf expression, with argument list

  /// cloning method
  expression *clone (Domain *d = NULL) const
    {return new exprOpp (argument_ -> clone (d));}

  /// the operator's function
  inline unary_function F () 
    {return opp;}

  void print (std::ostream &out, 
	      bool descend) const;

  /// print operator
  //std::string printOp () const
  //{return "-";}

  /// differentiation
  expression *differentiate (int index); 

  /// simplification
  virtual expression *simplify ();

  /// get a measure of "how linear" the expression is (see CouenneTypes.h)
  inline int Linearity ()
    {return argument_ -> Linearity ();}

  /// Get lower and upper bound of an expression (if any)
  void getBounds (expression *&, expression *&);

  /// special version for linear constraints
  virtual void generateCuts (expression *, const OsiSolverInterface &, 
			     OsiCuts &, const CouenneCutGenerator *,
			     t_chg_bounds * = NULL, int = -1, 
			     CouNumber = -COUENNE_INFINITY, 
			     CouNumber =  COUENNE_INFINITY);

  /// code for comparisons
  virtual enum expr_type code () 
    {return COU_EXPROPP;}

  /// is this expression integer?
  bool isInteger ()
    {return argument_ -> isInteger ();}

  /// implied bound processing
  bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *);

  /// standardization (to deal with complex arguments)
  exprAux *standardize (CouenneProblem *, bool addAux = true);
};

#endif

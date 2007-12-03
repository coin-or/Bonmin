/*
 * Name:    exprCopy.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the class exprCopy
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRCOPY_HPP
#define COUENNE_EXPRCOPY_HPP

#include <iostream>

#include "CouenneTypes.hpp"
#include "expression.hpp"

class CouenneObject;

// expression copy (points to VALUE of another expression) 

class exprCopy: public expression {

 protected:

  /// the expression this object is a (reference) copy of
  expression *copy_;

 public:

  /// node type
  inline enum nodeType Type () 
    {return copy_ -> Type ();}

  /// Constructor
  exprCopy  (expression *copy):
    copy_ (copy) {}

  /// Copy constructor
  exprCopy (const exprCopy &e) {
    copy_ = e.Original () -> clone ();
  }

  /// Cloning method
  virtual exprCopy *clone () const
    {return new exprCopy (*this);}

  /// If this is an exprClone of a exprClone of an expr???, point to
  /// the original expr??? instead of an exprClone -- improves computing
  /// efficiency
  inline const expression *Original () const
    {return copy_ -> Original ();}

  /// Get variable index in problem
  inline int Index () const
    {return copy_ -> Index ();}

  /// Return number of arguments (when applicable, that is, with N-ary functions)
  virtual inline int nArgs () const
    {return copy_ -> nArgs ();}

  /// return arglist (when applicable, that is, with N-ary functions)
  virtual inline expression **ArgList () const
    {return copy_ -> ArgList ();}

  /// return argument (when applicable, i.e., with univariate functions)
  virtual inline expression *Argument () const
    {return copy_ -> Argument ();}

  /// return pointer to argument (when applicable, i.e., with univariate functions)
  virtual inline expression **ArgPtr ()
    {return copy_ -> ArgPtr ();}

  /// I/O
  virtual void print (std::ostream &out = std::cout, 
		      bool descend      = false, 
		      CouenneProblem *p = NULL) const
    {copy_ -> Original () -> print (out, descend, p);}

  /// value (empty)
  virtual inline CouNumber Value () const 
    //    {return currValue_;}
    {return copy_ -> Value ();} // *** Check this! Should be the commented one 

  // TODO: FIX ME! a copy should just return an already evaluated
  // number, that's why it is very important that exprCopy should only
  // be used in successive evaluations.

  /// null function for evaluating the expression
  virtual inline CouNumber operator () () 
    {return (currValue_ = (*copy_) ());}
  //    {return (currValue_ = copy_ -> Value ());}

  /// differentiation
  inline expression *differentiate (int index) 
    {return copy_ -> differentiate (index);}

  /// fill in the set with all indices of variables appearing in the
  /// expression
  inline int DepList (std::set <int> &deplist, 
		      enum dig_type   type = ORIG_ONLY,
		      CouenneProblem *p    = NULL)
    {return copy_ -> DepList (deplist, type, p);}

  /// simplify expression (useful for derivatives)
  inline expression *simplify () 
    {return copy_ -> simplify ();}

  /// get a measure of "how linear" the expression is (see CouenneTypes.h)
  inline int Linearity ()
    {return copy_ -> Linearity ();}

  virtual inline bool isInteger ()
    {return copy_ -> isInteger ();}

  /// Get lower and upper bound of an expression (if any)
  inline void getBounds (expression *&lower, expression *&upper) 
    {copy_ -> getBounds (lower, upper);}

  /// Create standard formulation of this expression
  inline exprAux *standardize (CouenneProblem *p, bool addAux = true)
    {return copy_ -> standardize (p, addAux);}

  /// generate convexification cut for constraint w = this
  inline void generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg, 
			    t_chg_bounds *chg = NULL, int wind= -1, 
			    CouNumber lb = -COUENNE_INFINITY, 
			    CouNumber ub =  COUENNE_INFINITY)

    {copy_ -> generateCuts (w, si, cs, cg, chg, wind, lb, ub);}

  /// return an index to the variable's argument that is better fixed
  /// in a branching rule for solving a nonconvexity gap
  expression *getFixVar () 
    {return copy_ -> getFixVar ();}

  /// code for comparisons
  enum expr_type code () 
    {return copy_ -> code ();}

  /// either CONVEX, CONCAVE, AFFINE, or NONCONVEX 
  virtual enum convexity convexity ()
  {return copy_ -> convexity ();}

  /// compare this with other expression
  int compare (expression &e) 
    {return copy_ -> compare (e);}

  /// used in rank-based branching variable choice
  int rank (CouenneProblem *p)
    {return copy_ -> rank (p);} 

  /// implied bound processing
  bool impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg)
    {return copy_ -> impliedBound (wind, l, u, chg);}

  /// multiplicity of a variable: how many times this variable occurs
  /// in expressions throughout the problem
  virtual int Multiplicity ()
    {return copy_ -> Multiplicity ();}

  /// Set up branching object by evaluating many branching points for each expression's arguments.
  /// Return estimated improvement in objective function 
  virtual CouNumber selectBranch (const CouenneObject *obj,
				  const OsiBranchingInformation *info,
				  int     &ind,
				  double *&brpts,
				  int     &way) 
  {return copy_ -> selectBranch (obj, info, ind, brpts, way);}

  /// replace occurrence of a variable with another variable
  void replace (exprVar *, exprVar *);

  /// fill in dependence structure
  void fillDepSet (std::set <DepNode *, compNode> *dep, DepGraph *g)
    {copy_ -> fillDepSet (dep, g);}
};

#endif

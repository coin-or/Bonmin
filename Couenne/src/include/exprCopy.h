/*
 * Name:    exprCopy.h
 * Author:  Pietro Belotti
 * Purpose: definition of the class exprCopy
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRCOPY_H
#define COUENNE_EXPRCOPY_H

#include <iostream>

#include <CouenneTypes.h>
#include <expression.h>


// expression copy (points to VALUE of another expression) 

class exprCopy: public expression {

 protected:

  // the expression this object is a (reference) copy of
  expression *copy_;

 public:

  // node type
  inline enum nodeType Type () 
    {return copy_ -> Type ();}

  // Constructor, destructor
  exprCopy  (expression *copy):
    copy_ (copy) {}

  // copy constructor
  exprCopy (const exprCopy &e) {
    copy_ = e.Original () -> clone ();
  }

  // cloning method
  virtual exprCopy *clone () const
    {return new exprCopy (*this);}

  // If this is an exprClone of a exprClone of an expr???, point to
  // the original expr??? instead of an exprClone -- improves computing
  // efficiency
  inline const expression *Original () const
    {return copy_ -> Original ();}

  // get variable index in problem
  inline int Index () const
    {return copy_ -> Index ();}

  // I/O
  virtual void print (std::ostream &out) const
  {out << "["; copy_ -> Original () -> print (out); out << "]";}

  // value (empty)
  virtual inline CouNumber Value () const 
    //    {return currValue_;}
    {return copy_ -> Value ();} // *** Check this! Should be the commented one 

  // FIX ME! a copy should just return an already evaluated number,
  // that's why it is very important that exprCopy should only be used
  // in successive evaluations. 

  // null function for evaluating the expression
  virtual inline CouNumber operator () () 
    {return (currValue_ = (*copy_) ());}
  //    {return (currValue_ = copy_ -> Value ());}

  // differentiation
  inline expression *differentiate (int index) 
    {return copy_ -> differentiate (index);}

  // dependence on variable set
  inline bool dependsOn (int *varlist, int n) 
    {return copy_ -> dependsOn (varlist, n);}

  // simplify expression (useful for derivatives)
  inline expression *simplify () 
    {return copy_ -> simplify ();}

  // get a measure of "how linear" the expression is (see CouenneTypes.h)
  inline int Linearity ()
    {return copy_ -> Linearity ();}

  // Get lower and upper bound of an expression (if any)
  inline void getBounds (expression *&lower, expression *&upper) 
    {copy_ -> getBounds (lower, upper);}

  // Create standard formulation of this expression
  inline exprAux *standardize (CouenneProblem *p)
    {return copy_ -> standardize (p);}

  // generate convexification cut for constraint w = this
  inline void generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) 
    {copy_ -> generateCuts (w, si, cs, cg);}

  // return an index to the variable's argument that is better fixed
  // in a branching rule for solving a nonconvexity gap
  virtual expression *getFixVar () 
    {return copy_ -> getFixVar ();}

  ///
  virtual enum expr_type code () 
    {return copy_ -> code ();}

  ///
  int compare (expression &e) 
    {return copy_ -> compare (e);}

  /// used in rank-based branching variable choice
  virtual int rank (CouenneProblem *p)
    {return copy_ -> rank (p);} 

  /// implied bound processing
  bool impliedBound (int wind, CouNumber *l, CouNumber *u, char *chg)
    {return copy_ -> impliedBound (wind, l, u, chg);}
};

#endif

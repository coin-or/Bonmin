/*
 * Name:    exprVar.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the class exprVar for variables 
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRVAR_H
#define COUENNE_EXPRVAR_H

#include <iostream>
#include <set>

#include <CouenneTypes.h>
#include <expression.hpp>
#include <exprConst.hpp>
//#include <exprBound.hpp>

class CouenneProblem;


/// variable-type operator. All variables of the expression must be
/// objects of this class

class exprVar: public expression {

 protected:

  int varIndex_; //< the index of the variable's current value

 public:

  /// node type
  virtual inline enum nodeType Type () 
    {return VAR;}

  /// Constructor
  exprVar (int varIndex):
    varIndex_ (varIndex) {}

  /// destructor
  virtual ~exprVar () {}

  /// copy constructor
  exprVar (const exprVar &e):
    varIndex_ (e.Index ()) {}

  /// cloning method
  virtual exprVar *clone () const
    {return new exprVar (*this);}

  /// get variable index in problem
  inline int Index () const
    {return varIndex_;}

  /// for compatibility with exprAux
  virtual inline expression *Image () const
    {return NULL;}

  /// Bound get
  virtual expression *Lb ();
  virtual expression *Ub ();

  /// print
  virtual void print (std::ostream &out = std::cout, bool = false, CouenneProblem * = NULL) const
    {out << "x_" << varIndex_;}

  /// return the value of the variable
  virtual inline CouNumber operator () () 
    {return (currValue_ = expression::variables_ [varIndex_]);}

  /// return the value of the variable
  inline CouNumber Value ()
    {return currValue_;}

  /// differentiation
  virtual inline expression *differentiate (int index) 
    {return new exprConst ((index == varIndex_) ? 1 : 0);}

  /// dependence on variable set
  virtual int dependsOn (int *, int);

  /// set bounds depending on both branching rules and propagated
  /// bounds. To be used after standardization
  virtual inline void crossBounds () {}

  /// simplify
  inline expression *simplify () 
    {return NULL;}

  /// get a measure of "how linear" the expression is (see CouenneTypes.h)
  virtual inline int Linearity ()
    {return LINEAR;}

  /// is this expression integer?
  virtual bool isInteger ()
    {return false;}

  /// Get lower and upper bound of an expression (if any)
  virtual void getBounds (expression *&, expression *&);

  /// generate cuts for expression associated with this auxiliary
  virtual void generateCuts (const OsiSolverInterface &, 
			     OsiCuts &, const CouenneCutGenerator *, 
			     t_chg_bounds * = NULL, int = -1, 
			     CouNumber = -COUENNE_INFINITY, 
			     CouNumber =  COUENNE_INFINITY) {}

  /// general function to be specialized
  /*virtual void generateCuts (const OsiSolverInterface &, 
			     OsiCuts &, const CouenneCutGenerator *, 
			     t_chg_bounds *, int,
			     CouNumber, CouNumber) {}*/

  /// generate convexification cut for constraint w = this
  virtual void generateCuts (exprAux *w, const OsiSolverInterface &si, 
			     OsiCuts &cs, const CouenneCutGenerator *cg, 
			     t_chg_bounds * = NULL, int = -1, 
			     CouNumber = -COUENNE_INFINITY, 
			     CouNumber =  COUENNE_INFINITY);

  /// return an index to the variable's argument that is better fixed
  /// in a branching rule for solving a nonconvexity gap
  virtual expression *getFixVar () {return this;}

  /// code for comparison
  virtual enum expr_type code () {return COU_EXPRVAR;}

  /// implied bound processing
  virtual bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *);

  /// rank of an original variable is always one
  virtual int rank (CouenneProblem *p) 
    {return 1;}

  /// update dependence set with index of this variable
  virtual void fillDepSet (std::set <DepNode *, compNode> *, DepGraph *);
};

#endif

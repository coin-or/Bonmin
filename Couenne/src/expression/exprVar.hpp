/*
 * Name:    exprVar.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the class exprVar for variables 
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRVAR_HPP
#define COUENNE_EXPRVAR_HPP

#include <iostream>
#include <set>

#include "CouenneTypes.hpp"
#include "expression.hpp"
#include "exprConst.hpp"


/// variable-type operator
///
/// All variables of the expression must be objects of this class or
/// of the derived exprAux class

class exprVar: public expression {

 protected:

  int varIndex_; ///< the index of the variable's current value

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
  //virtual inline expression *Image () const
  //{return NULL;}

  // Bounds
  virtual expression *Lb (); ///< get lower bound expression
  virtual expression *Ub (); ///< get upper bound expression

  // Bounds
  virtual CouNumber lb () {return expression::lbounds_ [varIndex_];} ///< lower bound
  virtual CouNumber ub () {return expression::ubounds_ [varIndex_];} ///< upper bound

  /// print
  virtual void print (std::ostream &out = std::cout,
		      bool = false) const
    {out << "x_" << varIndex_;}

  /// return the value of the variable
  virtual inline CouNumber operator () () 
    {return expression::variables_ [varIndex_];}

  /// differentiation
  virtual inline expression *differentiate (int index) 
    {return new exprConst ((index == varIndex_) ? 1. : 0.);}

  /// fill in the set with all indices of variables appearing in the
  /// expression
  virtual inline int DepList (std::set <int> &deplist, 
			      enum dig_type type = ORIG_ONLY) {

    if (deplist.find (varIndex_) == deplist.end ()) {
      deplist.insert (varIndex_); 
      return 1;
    }
    return 0;
  }

  /// set bounds depending on both branching rules and propagated
  /// bounds. To be used after standardization
  virtual inline void crossBounds () {}

  /// simplify
  virtual inline expression *simplify () 
    {return NULL;}

  /// get a measure of "how linear" the expression is (see CouenneTypes.h)
  virtual inline int Linearity ()
    {return LINEAR;}

  /// is this variable integer?
  virtual inline bool isInteger () {
    CouNumber lb = (*(Lb ())) (); 
    return ((fabs (lb - (*(Ub ())) ()) < COUENNE_EPS) && (::isInteger (lb)));
  }

  /// Get lower and upper bound of an expression (if any)
  virtual void getBounds (expression *&, expression *&);

  /// generate cuts for expression associated with this auxiliary
  virtual void generateCuts (const OsiSolverInterface &, 
			     OsiCuts &, const CouenneCutGenerator *, 
			     t_chg_bounds * = NULL, int = -1, 
			     CouNumber = -COUENNE_INFINITY, 
			     CouNumber =  COUENNE_INFINITY) {}

  /// generate convexification cut for constraint w = this
  virtual void generateCuts (expression *w, 
			     const OsiSolverInterface &si, 
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
  virtual int rank () 
    {return 1;}

  /// update dependence set with index of this variable
  virtual void fillDepSet (std::set <DepNode *, compNode> *, DepGraph *);

  /// is this variable fixed?
  virtual inline bool isFixed ()
  {return (fabs ((*(Lb ())) () - (*(Ub ())) ()) < COUENNE_EPS);}
};

#endif

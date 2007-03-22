/*
 * Name:    exprAux.h
 * Author:  Pietro Belotti
 * Purpose: definition of the auxiliary variable class (used in
 *          standardization and convexification)
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRAUX_H
#define COUENNE_EXPRAUX_H

#include <iostream>
#include <exprVar.h>
#include <exprBound.h>
#include <exprMax.h>
#include <exprMin.h>
#include <CouenneTypes.h>
#include <CglCutGenerator.hpp>


// expression base class

class exprAux: public exprVar {

 protected:

  /// The expression associated with this auxiliary variable
  expression *image_;

  /// bounds determined by the associated expression's bounds and the
  /// relative constraint's bound
  expression *lb_;
  expression *ub_;

  /// used in rank-based branching variable choice: original variables
  /// have rank 1; auxiliary w=f(x) has rank r(w) = r(x)+1; finally,
  /// auxiliary w=f(x1,x2...,xk) has rank r(w) = 1+max{r(xi):i=1..k}.
  int rank_;

 public:

  // Node type
  inline enum nodeType Type () 
    {return AUX;}

  // Constructor
  exprAux (expression *, int, int);

  // Destructor
  ~exprAux () {
    delete image_; 
    delete lb_; 
    delete ub_;
  }

  // copy constructor
  exprAux (const exprAux &e):
    exprVar (e.varIndex_),
    image_  (e.Image () -> clone ()),
    rank_   (e.rank_)

    {image_ -> getBounds (lb_, ub_);}

  // cloning method
  virtual exprAux *clone () const
    {return new exprAux (*this);}

  // set lower bound
  void setLB (expression *lb) 
    {if (lb_) delete lb_; lb_ = lb;}

  // set upper bound
  void setUB (expression *ub) 
    {if (ub_) delete ub_; ub_ = ub;}

  // Bound get
  expression *Lb () {return lb_;}
  expression *Ub () {return ub_;}

  // I/O
  void print (std::ostream &out) const
    {out << "w_" << varIndex_;}

  // The expression associated with this auxiliary variable
  inline expression *Image () const
    {return image_;}

  // Null function for evaluating the expression
  inline CouNumber operator () () 
    {return (currValue_ = expression::Variable (varIndex_));}

  // Differentiation
  inline expression *differentiate (int index) 
    {return image_ -> differentiate (index);}

  // Dependence on variable set
  inline bool dependsOn (int *indices, int num) 
    {return image_ -> dependsOn (indices, num);}

  // Get a measure of "how linear" the expression is (see CouenneTypes.h)
  inline int Linearity ()
    {return LINEAR;
    /*return image_ -> Linearity ();*/}

  // Get lower and upper bound of an expression (if any)
  inline void getBounds (expression *&lb, expression *&ub) {
    expression *lb0, *ub0;
    image_ -> getBounds (lb0, ub0);
    lb = new exprMax (lb0, new exprLowerBound (varIndex_)); 
    ub = new exprMin (ub0, new exprUpperBound (varIndex_));
  }

  // generate cuts for expression associated with this auxiliary
  void generateCuts (const OsiSolverInterface &, 
		     OsiCuts &, const CouenneCutGenerator *);

  /// used in rank-based branching variable choice
  virtual int rank (CouenneProblem *p)
    {return rank_;} 
};

#endif

/*
 * Name:    CouenneProblemElem.hpp
 * Author:  Pietro Belotti
 * Purpose: define the classes used by class CouenneProblem
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_PROBLEM_ELEM_HPP
#define COUENNE_PROBLEM_ELEM_HPP

#include <iostream>

#include "CouenneTypes.hpp"
#include "expression.hpp"
#include "exprConst.hpp"


/** Class to represent nonlinear constraints
 *
 *  It consists of an expression as the body and two range expressions
 *  as lower- and upper bounds.  
 *
 *  A general constraint is defined as lb_ <= body_ <= ub_, where all
 *  three components are expressions, depending on variables,
 *  auxiliaries and bounds. If the constraint is 2 <= exp (x1+x2) <=
 *  4, then:
 *
 *  body_ = exp (x1+x2), that is, 
 *
 *  new exprExp (new exprSum (new exprVar (1), new exprVar (2))
 *
 *  while lb_ = new exprConst (2.) and ub_ = new exprConst (4.).
 */

class CouenneConstraint {

 protected:

  expression *body_; ///< body of constraint

  expression *lb_;   ///< lower bound
  expression *ub_;   ///< upper bound

 public:

  /// constructor
  CouenneConstraint  (expression *body = NULL, 
  	              expression *lb   = NULL, 
		      expression *ub   = NULL):
    body_     (body), 
    lb_       (lb), 
    ub_       (ub) {

    if (!lb_) 
      if (!ub_) {
	lb_ = new exprConst (0.);
	ub_ = new exprConst (0.);
      } 
      else         lb_ = new exprConst (- COUENNE_INFINITY);
    else if (!ub_) ub_ = new exprConst   (COUENNE_INFINITY);
  }

  /// destructor
  ~CouenneConstraint () {
    delete body_; 
    delete lb_; 
    delete ub_;
  }

  /// copy constructor
  CouenneConstraint  (const CouenneConstraint &c):
    body_  (c.Body () -> clone ()), 
    lb_    (c.Lb   () -> clone ()),
    ub_    (c.Ub   () -> clone ()) {}

  /// cloning method
  inline CouenneConstraint *clone () const
    {return new CouenneConstraint (*this);}

  // get constraint's elements
  inline expression *Lb   () const {return lb_;}   ///< expression of lower bound
  inline expression *Ub   () const {return ub_;}   ///< expression of upper bound
  inline expression *Body () const {return body_;} ///< expression of body of constraint

  /// set body of constraint
  inline expression *Body (expression *newBody) 
    {body_ = newBody; return body_;}

  /// decompose body of constraint through auxiliary variables
  exprAux *standardize (CouenneProblem *);

  /// print constraint
  void print (std::ostream & = std::cout);
};



/**
 * Objective function
 *
 * It consists of an expression and an optimization direction.
 */

class CouenneObjective {

 protected:

  /// expression to optimize
  expression *body_;

  /// can be MAXIMIZE or MINIMIZE
  enum opt_sense sense_;

 public:

  /// constructor
  CouenneObjective  (expression *body, enum opt_sense sense):
    body_ (body), sense_ (sense) {}

  /// destructor
  ~CouenneObjective () 
    {delete body_;}

  /// copy constructor
  CouenneObjective  (const CouenneObjective &o):
    body_  (o.body_ -> clone ()), 
    sense_ (o.sense_) {}

  /// cloning method
  inline CouenneObjective *clone () const
    {return new CouenneObjective (*this);}

  /// get optimization sense
  inline enum opt_sense Sense () const
    {return sense_;}

  /// get body
  inline expression *Body () const
    {return body_;}

  /// Set body
  expression *Body (expression *newBody) 
    {body_ = newBody; return body_;}

  /// Get standard form of this objective function
  inline exprAux *standardize (CouenneProblem *p) 
    {return body_ -> standardize (p);}

  /// Print to iostream
  void print (std::ostream &out = std::cout) {
    out << (sense_ == MAXIMIZE ? "max " : "min ");
    body_ -> print (out);
    out << std::endl;
  }
};

#endif

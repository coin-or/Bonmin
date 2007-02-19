/*
 * Name:    CouenneProblemElem.h
 * Author:  Pietro Belotti
 * Purpose: define the classes used by class CouenneProblem
 *
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_PROBLEM_ELEM_H
#define COUENNE_PROBLEM_ELEM_H

#include <iostream>

#include <CouenneTypes.h>
#include <expression.h>


// CouenneConstraint class, with any expression as the body and two range
// expressions

class CouenneConstraint {

 protected:

  // a general constraint is defined as lb_ <= body_ <= ub_, where all
  // three components are expressions, depending on variables,
  // auxiliaries and bounds. If the constraint is 2 <= exp (x1+x2) <=
  // 4, then:

  // body_ = exp (x1+x2), that is, 
  // new exprExp (new exprSum (new exprVar (1), new exprVar (2)),
  expression *body_;

  // while lb_ = new exprConst (2) and ub_ = new exprConst (4).
  expression *lb_;
  expression *ub_;

 public:

  // constructor
  CouenneConstraint  (expression *body = NULL, 
  	              expression *lb   = NULL, 
		      expression *ub   = NULL):
    body_ (body), 
    lb_   (lb), 
    ub_   (ub) {

    if (!lb_) 
      if (!ub_) {
	lb_ = new exprConst (0);
	ub_ = new exprConst (0);
      } 
      else         lb_ = new exprConst (- COUENNE_INFINITY);
    else if (!ub_) ub_ = new exprConst   (COUENNE_INFINITY);
  }

  // destructor
  ~CouenneConstraint () {
    delete body_; 
    delete lb_; 
    delete ub_;
  }

  // copy constructor
  CouenneConstraint  (const CouenneConstraint &c):
    body_  (c.Body () -> clone ()), 
    lb_    (c.Lb   () -> clone ()),
    ub_    (c.Ub   () -> clone ()) {}

  // cloning method
  inline CouenneConstraint *clone () const
    {return new CouenneConstraint (*this);}

  // get constraint's elements
  inline expression *Lb   () const {return lb_;}
  inline expression *Ub   () const {return ub_;}
  inline expression *Body () const {return body_;}

  // set body of constraint
  inline expression *Body (expression *newBody) 
    {body_ = newBody; return body_;}

  // decompose body of constraint through auxiliary variables
  inline exprAux *standardize (CouenneProblem *p) 
    {return body_ -> standardize (p);}

  // print constraint
  void print (std::ostream &);
};


// Objective function class, with an expression and an optimization
// direction

class Objective {

 protected:

  // expression to optimize
  expression *body_;

  // can be COUENNE_MAXIMIZE or COUENNE_MINIMIZE
  enum opt_sense sense_;

 public:

  // constructor
  Objective  (expression *body, enum opt_sense sense):
    body_ (body), sense_ (sense) {}

  // destructor
  ~Objective () 
    {delete body_;}

  // copy constructor
  Objective  (const Objective &o):
    body_  (o.Body  () -> clone ()), 
    sense_ (o.Sense ()) {}

  // cloning method
  inline Objective *clone () const
    {return new Objective (*this);}

  // optimization sense
  inline enum opt_sense Sense () const
    {return sense_;}

  // get body
  inline expression *Body () const
    {return body_;}

  // set body
  expression *Body (expression *newBody) 
    {body_ = newBody; return body_;}

  // get standard form of this objective function
  inline exprAux *standardize (CouenneProblem *p) 
    {return body_ -> standardize (p);}

  // I/O
  void print (std::ostream &out = std::cout) {
    out << (sense_ == MAXIMIZE ? "max " : "min ");
    body_ -> print (out);
    out << std::endl;
  }
};

#endif

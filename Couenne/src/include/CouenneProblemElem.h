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
  // three components are general expressions, depending on variables,
  // auxiliaries and bounds. If the constraint is 2 <= exp (x1+x2) <=
  // 4, then:

  // body_ = exp (x1+x2), that is, new exprExp (new exprSum (new
  // exprVar (1), new exprVar (2)),
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
    ub_   (ub) 
    {}

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


// Linear constraint class, with an array of expression for the
// coefficients, an array of indices (int) for the variables, a lower
// bound and an upper bound (both of class expression). 
/*
class LinearConstraint {

 protected:

  // number of coefficients
  int nterms_;

  // array of coefficients (all of them are expressions)
  expression **coeff_;

  // array of indices of variables/auxiliaries associated with coeff
  int *indices_;

  // right-hand side
  expression *rhs_;

  // can be COUENNE_EQ, COUENNE_RNG, COUENNE_GE, COUENNE_LE
  enum con_sign sign_;

 public:

  int            nTerms  () {return nterms_;}
  expression   **Coeff   () {return coeff_;}
  expression    *Rhs     () {return rhs_;}
  int           *Indices () {return indices_;}
  enum con_sign  Sign    () {return sign_;}

  // general constructor: specify array of coefficients and indices,
  // lower- and upper bound, and sign
  LinearConstraint  (int nterms         = 0,
		     expression **coeff = NULL, 
		     int *indices       = NULL,
		     expression *rhs    = NULL,
		     enum con_sign sign = COUENNE_EQ):

    nterms_  (nterms),
    coeff_   (coeff), 
    indices_ (indices),
    rhs_     (rhs),
    sign_    (sign) 
    {
      if (!rhs_) 
	rhs_ = new exprConst (0);
    }

  // constructor specific for auxiliary variables: specify auxiliary
  // variable w, argument of unary function x, its coefficients and
  // right hand side in the constraint -w + ax >=< b
  LinearConstraint  (exprAux *w,
		     exprVar *x,
		     expression *coeff,
		     expression *rhs = NULL,
		     enum con_sign sign = COUENNE_EQ):

    nterms_  (2),
    coeff_   (new expression * [2]), 
    indices_ (new int [2]),
    rhs_     (rhs), 
    sign_    (sign) 
    {
      if (!rhs_) 
	rhs_ = new exprConst (0);

      coeff_ [0] = new exprConst (-1);    indices_ [0] = w -> Index ();
      coeff_ [1] = coeff;                 indices_ [1] = x -> Index ();

      coeff_ [0] -> print (std::cout); printf ("  ");
      coeff_ [1] -> print (std::cout); printf ("\n");
    }

  // destructor
  ~LinearConstraint () {
    while (nterms_--) 
      delete coeff_ [nterms_];
    delete [] coeff_;
    delete [] indices_;
    if (rhs_) delete rhs_;
  }

  void print (std::ostream &);
};
*/

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

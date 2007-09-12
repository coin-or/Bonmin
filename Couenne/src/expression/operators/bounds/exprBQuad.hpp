/*
 * Name:    exprBQuad.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of operators to compute lower/upper bounds of quadratic forms
 *
 * (C) Pietro Belotti 2006. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRBQUAD_H
#define COUENNE_EXPRBQUAD_H

#include <exprOp.hpp>
#include <exprQuad.hpp>


/// method to actually compute the bound
CouNumber computeQBound (int, exprQuad *);


/// class to compute lower bound of a fraction based on the bounds of
/// both numerator and denominator

class exprLBQuad: public expression {

  exprQuad *ref_; //< quadratic form, reference expression

 public:

  /// Constructor
  exprLBQuad (exprQuad *ref): ref_ (ref) {}

  /// copy constructor
  exprLBQuad (const exprLBQuad &src): 
    ref_ (dynamic_cast <exprQuad *> (src.ref_ -> clone ())) {}

  /// destructor
  ~exprLBQuad () {}

  /// cloning method
  expression *clone () const
    {return new exprLBQuad (*this);}

  /// function for the evaluation of the expression
  CouNumber operator () () 
  {return computeQBound (-1, NULL);}

  /// I/O
  virtual void print (std::ostream &s = std::cout,     //< output stream
		      bool descend = false,            //< descend into auxiliaries' image?
		      CouenneProblem *p = NULL) const  //< problem pointer (in exprGroup)

  {s << "quadLower("; ref_ -> print (s, descend, p); s << ')';}
};


/// class to compute upper bound of a fraction based on the bounds of
/// both numerator and denominator

class exprUBQuad: public expression {

  exprQuad *ref_; //< quadratic form, reference expression

 public:

  /// Constructor
  exprUBQuad (exprQuad *ref): ref_ (ref) {}

  /// copy constructor
  exprUBQuad (const exprUBQuad &src): 
    ref_ (dynamic_cast <exprQuad *> (src.ref_ -> clone ())) {}

  /// destructor
  ~exprUBQuad () {}

  /// cloning method
  expression *clone () const
    {return new exprUBQuad (*this);}

  /// function for the evaluation of the expression
  CouNumber operator () () 
  {return computeQBound (1, NULL);}

  /// I/O
  virtual void print (std::ostream &s = std::cout,     //< output stream
		      bool descend = false,            //< descend into auxiliaries' image?
		      CouenneProblem *p = NULL) const  //< problem pointer (in exprGroup)

  {s << "quadUpper("; ref_ -> print (s, descend, p); s << ')';}
};

#endif

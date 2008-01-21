/*
 * Name:    exprBQuad.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of operators to compute lower/upper bounds of quadratic forms
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRBQUAD_H
#define COUENNE_EXPRBQUAD_H

#include "exprOp.hpp"
#include "exprQuad.hpp"


/// class to compute lower bound of a fraction based on the bounds of
/// both numerator and denominator

class exprLBQuad: public expression {

  exprQuad *ref_; ///< quadratic form, reference expression

 public:

  /// Constructor
  exprLBQuad (exprQuad *ref): 
    ref_ (ref) {}

  /// copy constructor
  exprLBQuad (const exprLBQuad &src, const std::vector <exprVar *> *variables = NULL): 
    ref_ (dynamic_cast <exprQuad *> (src.ref_ -> clone (variables))) {}

  /// destructor
  ~exprLBQuad () {}

  /// cloning method
  expression *clone (const std::vector <exprVar *> *variables = NULL) const
    {return new exprLBQuad (*this, variables);}

  /// function for the evaluation of the expression
  inline CouNumber operator () () 
  {return ref_ -> computeQBound (-1);}

  /// I/O
  virtual void print (std::ostream &s = std::cout,     //< output stream
		      bool descend = false) const      //< descend into auxiliaries' image?

  {s << "quadLower("; ref_ -> print (s, descend); s << ')';}
};


/// class to compute upper bound of a fraction based on the bounds of
/// both numerator and denominator

class exprUBQuad: public expression {

  exprQuad *ref_; ///< quadratic form, reference expression

 public:

  /// Constructor
  exprUBQuad (exprQuad *ref): 
    ref_ (ref) {}

  /// copy constructor
  exprUBQuad (const exprUBQuad &src, const std::vector <exprVar *> *variables = NULL): 
    ref_ (dynamic_cast <exprQuad *> (src.ref_ -> clone (variables))) {}

  /// destructor
  ~exprUBQuad () {}

  /// cloning method
  expression *clone (const std::vector <exprVar *> *variables = NULL) const
    {return new exprUBQuad (*this, variables);}

  /// function for the evaluation of the expression
  inline CouNumber operator () () 
  {return ref_ -> computeQBound (1);}

  /// I/O
  virtual void print (std::ostream &s = std::cout,     //< output stream
		      bool descend = false) const      //< descend into auxiliaries' image?

  {s << "quadUpper("; ref_ -> print (s, descend); s << ')';}
};

#endif

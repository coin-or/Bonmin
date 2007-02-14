/*
 * Name:    exprOp.h
 * Author:  Pietro Belotti
 * Purpose: definition of the n-ary expression class
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPROP_H
#define COUENNE_EXPROP_H

#include <iostream>

#include <expression.h>
#include <CouenneTypes.h>


// general n-ary operator-type expression: requires argument list. All
// non-unary and non-leaf operators, i.e., sum, subtraction,
// multiplication, power, division, max, min, etc. are derived from
// this class.
//

class exprOp: public expression {

 protected:

  expression **arglist_; //< argument list is an array of pointers to other expressions
  int          nargs_;   //< number of arguments (cardinality of arglist)

 public:

  // node type
  virtual inline enum nodeType Type () 
    {return N_ARY;}

  // Constructors, destructor
  exprOp (expression **arglist, int nargs):  //< non-leaf expression, with argument list 
    arglist_ (arglist),
    nargs_   (nargs)
    {}

  exprOp (expression *arg0, expression *arg1):  //< two arguments 
    arglist_ (new expression * [2]),
    nargs_   (2)
    {arglist_ [0] = arg0; arglist_ [1] = arg1;}

 ~exprOp ();

  // copy constructor
  exprOp (const exprOp &e):
    arglist_ (new expression * [e.nArgs ()]),
    nargs_   (e.nArgs ()) {

  }

  // cloning method
  virtual exprOp *clone ()
    {return new exprOp (*this);}


  // return class members
  inline expression **ArgList () const 
    {return arglist_;}

  inline int nArgs () const 
    {return nargs_;}

  // string equivalent
  virtual const std::string name () const;

  // I/O
  virtual void print (std::ostream &, const std::string &, enum pos) const;

  // function for the evaluation of the expression
  virtual inline CouNumber operator () ();

  // dependence on variable set
  virtual bool dependsOn (int *, int);

  // simplification
  virtual expression *simplify ();

  // clone argument list (for use with clone method
  expression **clonearglist () const {
    expression **al = new expression * [nargs_];
    for (int i=0; i<nargs_; i++)
      al [i] = arglist_ [i] -> clone ();
    return al;
  }

  // compress argument list
  int shrink_arglist (CouNumber, CouNumber);

  // get a measure of "how linear" the expression is:
  //
  // CONSTANT  = 0: a constant
  // LINEAR    = 1: linear
  // QUADRATIC = 2: quadratic
  // NONLINER  = 3: nonlinear non-quadratic
  virtual inline int Linearity ()
    {return NONLINEAR;}

  exprAux *standardize (CouenneProblem *);

  // return an index to the variable's argument that is better fixed
  // in a branching rule for solving a nonconvexity gap
  virtual int getFixIndex () {return arglist_ [0] -> Index ();}
};


// expression evaluation -- n-ary operator (non-variable, non-constant)

inline CouNumber exprOp::operator () () {

  // Fetch argument list and compute it "recursively" (the operator()
  // of the elements in the list is called) to fill in the vector
  // containing the numerical value of the argument list.

  register expression **al = arglist_;

  for (register int i = nargs_; i--;) 
    *++sp = (**al++) ();

  return 0;
}

#endif

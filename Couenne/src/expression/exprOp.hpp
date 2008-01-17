/*
 * Name:    exprOp.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of the n-ary expression class
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPROP_HPP
#define COUENNE_EXPROP_HPP

#include <iostream>

#include "expression.hpp"
#include "CouenneTypes.hpp"

class CouenneProblem;

/// general n-ary operator-type expression: requires argument
/// list. All non-unary and non-leaf operators, i.e., sum,
/// subtraction, multiplication, power, division, max, min, etc. are
/// derived from this class.

class exprOp: public expression {

 protected:

  expression **arglist_; ///< argument list is an array of pointers to other expressions
  int          nargs_;   ///< number of arguments (cardinality of arglist)

 public:

  /// Node type
  virtual inline enum nodeType Type () 
    {return N_ARY;}

  /// Constructor
  exprOp (expression **arglist, int nargs):  //< non-leaf expression, with argument list 
    arglist_ (arglist),
    nargs_   (nargs)
    {}

  /// Constructor with two arguments (for convenience)
  exprOp (expression *arg0, expression *arg1):  //< two arguments 
    arglist_ (new expression * [2]),
    nargs_   (2)
    {arglist_ [0] = arg0; arglist_ [1] = arg1;}

  /// Destructor
  virtual ~exprOp ();

  /// Copy constructor: only allocate space for argument list, which
  /// will be copied with clonearglist()
  exprOp (const exprOp &e):
    arglist_ (new expression * [e.nArgs ()]),
    nargs_   (e.nArgs ()) {}

  /// cloning method
  virtual expression *clone () const
    {return new exprOp (*this);}

  /// return argument list
  inline expression **ArgList () const 
    {return arglist_;}

  /// return number of arguments
  inline int nArgs () const 
    {return nargs_;}

  /// I/O
  virtual void print (std::ostream &out = std::cout,
		      bool = false) const;

  /// print position (PRE, INSIDE, POST)
  virtual enum pos printPos () const
    {return INSIDE;}

  /// print operator
  virtual std::string printOp () const
    {return "??";}

  /// fill in the set with all indices of variables appearing in the
  /// expression
  virtual inline int DepList (std::set <int> &deplist, 
			      enum dig_type type = ORIG_ONLY) {
    int tot = 0;
    for (int i = nargs_; i--;)
      tot += arglist_ [i] -> DepList (deplist, type);
    return tot;
  }

  /// simplification
  virtual expression *simplify ();

  /// clone argument list (for use with clone method
  expression **clonearglist () const {
    if (nargs_) {
      expression **al = new expression * [nargs_];
      for (register int i=0; i<nargs_; i++)
	al [i] = arglist_ [i] -> clone ();
      return al;
    } else return NULL;
  }

  /// compress argument list
  int shrink_arglist (CouNumber, CouNumber);

  /// get a measure of "how linear" the expression is (see CouenneTypes.h)
  virtual inline int Linearity ()
    {return NONLINEAR;}

  /// generate auxiliary variable
  virtual exprAux *standardize (CouenneProblem *, bool addAux = true);

  /// return an index to the variable's argument that is better fixed
  /// in a branching rule for solving a nonconvexity gap
  virtual expression *getFixVar ()
    {printf ("### Warning: called empty exprOp::getFixVar()\n"); return arglist_ [0];}

  /// return code to classify type of expression
  virtual enum expr_type code ()
    {return COU_EXPROP;}

  /// is this expression integer?
  virtual bool isInteger ();

  /// compare with other generic exprOp
  virtual int compare (exprOp &);

  /// used in rank-based branching variable choice
  virtual int rank ();

  /// fill in dependence structure
  /// update dependence set with index of this variable
  virtual void fillDepSet (std::set <DepNode *, compNode> *dep, DepGraph *g) {
    for (int i=nargs_; i--;)
      arglist_ [i] -> fillDepSet (dep, g);
  }

  /// replace variable with other
  virtual void replace (exprVar *, exprVar *);
};

#endif

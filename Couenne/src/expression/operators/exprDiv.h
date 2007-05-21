/*
 * Name:    exprDiv.h
 * Author:  Pietro Belotti
 * Purpose: definition of divisions
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRDIV_H
#define COUENNE_EXPRDIV_H

#include <exprOp.h>
#include <CouennePrecisions.h>


// class for divisions

class exprDiv: public exprOp {

 public:

  // Constructors, destructor
  exprDiv (expression **al, int n = 2): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  exprDiv (expression *arg0, expression *arg1):
    exprOp (arg0, arg1) {}

  // cloning method
  expression *clone () const
    {return new exprDiv (clonearglist (), nargs_);}

  /// print operator
  std::string printOp () const
    {return "/";}

  /// print position
  //  virtual enum pos printPos () 
  //    {return INSIDE;}

  // I/O
  //  void print (std::ostream&) const;

  // function for the evaluation of the expression
  inline CouNumber operator () ();

  // differentiation
  expression *differentiate (int index); 

  // simplification
  expression *simplify ();

  // get a measure of "how linear" the expression is (see CouenneTypes.h)
  inline int Linearity () {

    if (arglist_ [1] -> Type () == CONST)
         return arglist_ [0] -> Linearity ();
    else return NONLINEAR;
  }

  // Get lower and upper bound of an expression (if any)
  void getBounds (expression *&lb, expression *&ub);

  // reduce expression in standard form, creating additional aux
  // variables (and constraints)
  exprAux *standardize (CouenneProblem *p);

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);

  // return an index to the variable's argument that is better fixed
  // in a branching rule for solving a nonconvexity gap
  expression *getFixVar ();

  ///
  virtual enum expr_type code () {return COU_EXPRDIV;}

  /// implied bound processing
  bool impliedBound (int, CouNumber *, CouNumber *, char *);

  /// set up branching object by evaluating many branching points for
  /// each expression's arguments
  CouNumber selectBranch (expression *, const OsiBranchingInformation *,
			  int &, double * &, int &);

  /*  /// distance covered by current point if branching rule applied to this expression
  double BranchGain (expression *, const OsiBranchingInformation *);

  /// branching object best suited for this expression
  OsiBranchingObject *BranchObject (expression *, const OsiBranchingInformation *);*/
};


// compute division

inline CouNumber exprDiv::operator () () {

  //  exprOp:: operator () ();

  //  register CouNumber denominator = *sp--;
  //  return (currValue_ = (*sp-- / denominator));
  return (currValue_ = (*(arglist_ [0])) () / (*(arglist_ [1])) ());
}


#define SAFE_COEFFICIENT 1e9

// check if bounding box is suitable for a multiplication/division
// convexification constraint

inline bool is_boundbox_regular (register CouNumber b1, register CouNumber b2) {

  // Why SAFE_COEFFICIENT and not COUENNE_INFINITY? Because OsiRowCut::set[LU]b do
  // not work for values more than SAFE_COEFFICIENT and apparently makes the
  // convexification infeasible.
  return 
    (fabs (b1)    < SAFE_COEFFICIENT) && 
    (fabs (b2)    < SAFE_COEFFICIENT) && 
    (fabs (b1*b2) < SAFE_COEFFICIENT);
    //    && ((fabs (b1) > COUENNE_EPS) || (fabs (b2) > COUENNE_EPS));
}

#endif

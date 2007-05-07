/*
 * Name:    exprMax.h
 * Author:  Pietro Belotti
 * Purpose: definition of $\f(x_{\argmax_{i\in I} y_i})$ 
 *
 * (C) Pietro Belotti 2006. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRMAX_H
#define COUENNE_EXPRMAX_H

#include <exprOp.h>
#include <exprClone.h>


//  class max

class exprMax: public exprOp {

 public:

  // Constructors, destructor
  exprMax  (expression **al, int n): 
    exprOp (al, n) {} //< non-leaf expression, with argument list

  exprMax  (expression *el0, expression *el1):
    exprOp (new expression * [4], 4) {
    arglist_ [0] = el0; arglist_ [1] = new exprClone (el0);
    arglist_ [2] = el1; arglist_ [3] = new exprClone (el1);
  }

  // cloning method
  exprMax *clone () const
    {return new exprMax (clonearglist (), nargs_);}

  // I/O
  void print (std::ostream &out) const
    {exprOp:: print (out, "max", PRE);}

  // function for the evaluation of the expression
  inline CouNumber operator () ();

  // differentiation
  inline expression *differentiate (int) 
    {return NULL;} 

  // simplification
  inline expression *simplify () 
    {return NULL;}

  // get a measure of "how linear" the expression is (see CouenneTypes.h)
  virtual inline int Linearity () 
    {return NONLINEAR;}

  // Get lower and upper bound of an expression (if any)
  //  void getBounds (expression *&, expression *&);

  // reduce expression in standard form, creating additional aux
  // variables (and constraints)
  virtual inline exprAux *standardize (CouenneProblem *)
    {return NULL;}

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);

  ///
  virtual enum expr_type code () {return COU_EXPRMAX;}
};


// compute maximum

inline CouNumber exprMax::operator () () {

  CouNumber best_val = (*(arglist_ [0])) ();
  int best_ind = 0;

  for (int ind = 2; ind < nargs_; ind += 2) {

    register CouNumber val = (*(arglist_ [ind])) ();

    if (val > best_val) {
      best_ind = ind;
      best_val = val;
    }
  }

  best_val = (*(arglist_ [best_ind + 1])) ();

  return (currValue_ = best_val);

  /*
  exprOp:: operator () ();

  register CouNumber best_val = *sp--; 
  register CouNumber best_el  = *sp--; 

  int n = nargs_ / 2;

  while (--n) {

    printf ("::max: curr val = %.3f\n", *sp);

    register CouNumber val = *sp--;
    register CouNumber el  = *sp--;

    if (el > best_el) {
      best_el  = el;
      best_val = val;
    }
  }

  return (currValue_ = best_val);
  */
}

#endif

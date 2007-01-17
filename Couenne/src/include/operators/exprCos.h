/*
 * Name:    exprCos.h
 * Author:  Pietro Belotti
 * Purpose: definition of cosine 
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRCOS_H
#define COUENNE_EXPRCOS_H

#include <exprUnary.h>
#include <exprConst.h>
#include <CouennePrecisions.h>

#include <math.h>


// class cosine

class exprCos: public exprUnary {

 public:

  // constructor, destructor
  exprCos  (expression *al): 
    exprUnary (al, cos) {}

  ~exprCos () {}

  // cloning method
  expression *clone () const
    {return new exprCos (argument_ -> clone ());}

  // print "cos" and argument
  void print (std::ostream&);

  // obtain derivative of expression
  expression *differentiate (int index); 

  // Get lower and upper bound of an expression (if any)
  virtual inline void getBounds (expression *&lb, expression *&ub)
    {lb = new exprConst (-1); ub = new exprConst (1);}

  // construct linear under-estimator for expression within problem *p
  // (p is used to add convexification constraints)
  //  int lowerLinearHull (exprAux *, int *&, expression ***&, 
  //		       int **&, expression **&, enum con_sign *&);

  // construct linear over-estimator for expression within problem *p
  // (p is used to add convexification constraints)
  //  int upperLinearHull (exprAux *, int *&, expression ***&, 
  //		       int **&, expression **&, enum con_sign *&);

  // generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);
};


// common convexification method used by both cos and sin

void trigGenCuts (exprAux *, OsiCuts &, const CouenneCutGenerator *, unary_function);


// convex envelope for sine/cosine 

void addHexagon (const CouenneCutGenerator *, // pointer to the caller cut generator 
		 OsiCuts &,      // cut set to be enriched
		 unary_function, // sine or cosine
		 bool,           // should violation be checked before adding cut?
		 exprAux *,      // auxiliary variable
		 expression *);  // argument of cos/sin (should be a variable)

#endif

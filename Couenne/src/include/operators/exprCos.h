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


/// class cosine

class exprCos: public exprUnary {

 public:

  /// constructor, destructor
  exprCos  (expression *al): 
    exprUnary (al) {}

  /// cloning method
  expression *clone () const
    {return new exprCos (argument_ -> clone ());}

  //// the operator's function
  inline unary_function F () {return cos;}

  /// print "cos" and argument
  void print (std::ostream&) const;

  /// obtain derivative of expression
  expression *differentiate (int index); 

  /// Get lower and upper bound of an expression (if any)
  virtual inline void getBounds (expression *&lb, expression *&ub)
    {lb = new exprConst (-1); ub = new exprConst (1);}

  /// generate equality between *this and *w
  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg);

  /// code for comparisons
  virtual enum expr_type code () {return COU_EXPRCOS;}
};


/// common convexification method used by both cos and sin

//void trigGenCuts (exprAux *, OsiCuts &, const CouenneCutGenerator *, unary_function);


CouNumber trigNewton (CouNumber, int);

/// convex envelope for sine/cosine 

void addHexagon (const CouenneCutGenerator *, // pointer to the caller cut generator 
		 OsiCuts &,      // cut set to be enriched
		 unary_function, // sine or cosine
		 exprAux *,      // auxiliary variable
		 expression *);  // argument of cos/sin (should be a variable)

#endif

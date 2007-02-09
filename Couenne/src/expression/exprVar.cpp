/*
 * Name:    expression.C
 * Author:  Pietro Belotti
 * Purpose: methods of the expression class
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <expression.h>
#include <exprAux.h>
#include <exprOp.h>
#include <exprUnary.h>
#include <exprVar.h>
#include <exprBound.h>

#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>


// is the variable one of those in varlist?

bool exprVar::dependsOn (int *varlist = NULL, int n = 1) {

  if (!varlist) 
    return true;

  while (n--) 
    if (varIndex_ == *varlist++) 
      return true;

  return false;
}


// Get lower and upper bound of a variable expression (if any)

void exprVar::getBounds (expression *&lb, expression *&ub) {

  lb = new exprLowerBound (varIndex_); 
  ub = new exprUpperBound (varIndex_);
}


// generate convexification cut for constraint w = this

void exprVar::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {
  OsiRowCut *cut;

  if ((cut = cg -> createCut (0., 0, w -> Index (), 1., 
			      varIndex_, -1., -1, 0., true)))
    cs.insert (cut);
  /*
 = new OsiRowCut;

  CouNumber *coeff = new CouNumber [2];  
  int       *index = new int       [2];  

  cut -> setLb (0);
  cut -> setUb (0);

  coeff [0] =  1; index [0] = w -> Index ();
  coeff [1] = -1; index [1] = varIndex_;

  cs.insert (cut);
  */
}


// auxiliary expression Constructor

exprAux::exprAux (expression *image, int index): 

  exprVar (index),
  image_  (image) {

  image_ -> getBounds (lb_, ub_);
  //  lb_ = new exprConst (- COUENNE_INFINITY);
  //  ub_ = new exprConst (  COUENNE_INFINITY);
}

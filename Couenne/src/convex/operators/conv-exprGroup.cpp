/*
 * Name:    conv-exprGroup.cpp
 * Author:  Pietro Belotti
 * Purpose: implementation of convexification methods for exprGroup
 *
 * (C) Pietro Belotti 2007. This file is licensed under the Common Public License (CPL)
 */

#include <OsiRowCut.hpp>
#include <OsiCuts.hpp>

#include <exprGroup.h>
#include <exprBound.h>
#include <exprMul.h>

#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>

/// Get lower and upper bound of an expression (if any)
void exprGroup::getBounds (expression *&lb, expression *&ub) {

  expression *lbnl, *ubnl;

  // compute lower/upper bound of nonlinear part
  exprSum::getBounds (lbnl, ubnl);

  // count linear and constant terms
  int nlin = 0;
  if (fabs (c0_) > COUENNE_EPS) nlin++;
  for (register int *ind = index_; *ind++>=0; nlin++);

  expression 
    **linall = new expression * [nlin + 1], // linear arglist for lower bound
    **linalu = new expression * [nlin + 1]; //                    upper

  // add constant to bounds
  if (fabs (c0_) > COUENNE_EPS) {
    *linall++ = new exprConst (c0_);
    *linalu++ = new exprConst (c0_);
  }

  for (register int *ind = index_, i=0; *ind>=0;) {

    CouNumber coeff = coeff_ [i++];

    expression *l = new exprLowerBound (*ind),
               *u = new exprUpperBound (*ind++);

    if (fabs (coeff - 1.) < COUENNE_EPS) {
      *linall++ = l;
      *linalu++ = u;
    } else {

      expression *c1 = new exprConst (coeff),
                 *c2 = new exprConst (coeff);

      if (coeff < 0) {
	*linall++ = new exprMul (c1, u);
	*linalu++ = new exprMul (c2, l);
      } else {
	*linall++ = new exprMul (c1, l);
	*linalu++ = new exprMul (c2, u);
      }
    }
  }

  *linall = lbnl;
  *linalu = ubnl;

  lb = new exprSum (linall - nlin, nlin + 1);
  ub = new exprSum (linalu - nlin, nlin + 1);
}


/// reduce expression in standard form, creating additional aux
/// variables (and constraints)
exprAux *exprGroup::standardize (CouenneProblem *p) {

  // same as in exprSum... the only difference is that we also have to
  // arrange it with the linear and constant term later in generateCuts
  exprOp::standardize (p);

  // create auxiliary pointing to this expression
  return p -> addAuxiliary (this);
}


// generate equality between *this and *w
void exprGroup::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			      OsiCuts &cs, const CouenneCutGenerator *cg) {

  // very similar to exprSum::generateCuts. First of all, this has
  // been standardized into a sum, so it only gets a cut in the
  // initial relaxation
  if (!(cg -> isFirst ()))
    return;

  // second, count the number of linear terms
  int nlin = 0;

  for (register int *ind = index_; *ind++>=0; nlin++);

  CouNumber *coeff = new CouNumber [nargs_ + nlin + 1];
  int       *index = new int       [nargs_ + nlin + 1];
  OsiRowCut *cut   = new OsiRowCut;

  CouNumber rhs = c0_;

  // first, make room for aux variable
  coeff [0] = -1.; index [0] = w -> Index ();

  // now add linear terms
  for (register int i=0; i<nlin; i++) {

    coeff [i+1] = coeff_ [i]; 
    index [i+1] = index_ [i];
  }

  int nv = nlin+1;

  // scan arglist for (aux) variables and constants
  for (register int i=0; i<nargs_; i++) {

    expression *curr = arglist_ [i];

    if (curr -> Type () == CONST)
      rhs += curr -> Value ();
    else {
      coeff [nv]   = 1.; 
      index [nv++] = curr -> Index ();
    }
  }

  cut -> setRow (nv, index, coeff);

  rhs = - rhs;

  cut -> setUb (rhs);
  cut -> setLb (rhs);

  // added only once, it is global
  cut -> setGloballyValid ();

  cs.insert (cut);

  delete [] index;
  delete [] coeff;

  delete cut;
}

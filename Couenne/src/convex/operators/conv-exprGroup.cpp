/*
 * Name:    conv-exprGroup.cpp
 * Author:  Pietro Belotti
 * Purpose: implementation of convexification methods for exprGroup
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <OsiRowCut.hpp>
#include <OsiCuts.hpp>

#include <exprGroup.hpp>
#include <exprBound.hpp>
#include <exprMul.hpp>

#include <CouenneProblem.hpp>
#include <CouenneCutGenerator.hpp>

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


// generate equality between *this and *w
void exprGroup::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			      OsiCuts &cs, const CouenneCutGenerator *cg,
			      t_chg_bounds *chg, 
			      int wind, CouNumber lb, CouNumber ub) {

  // very similar to exprSum::generateCuts. First of all, this has
  // been standardized into a sum, so it only gets a cut in the
  // initial relaxation
  if (!(cg -> isFirst ()))
    return;

  // there is one linear term so far: -w
  int nterms = 0;

  OsiRowCut *cut = new OsiRowCut;

  // count terms in linear part
  for (register int *ind = index_; *ind++ >= 0; nterms++);

  int displacement = (wind < 0) ? 1: 0;

  CouNumber *coeff = new CouNumber [nargs_ + nterms + displacement];
  int       *index = new int       [nargs_ + nterms + displacement];

  if (wind < 0) {
    // first, make room for aux variable
    coeff [0] = -1.; 
    index [0] = w -> Index ();
    lb = ub = 0;
  }

  lb -= c0_;
  ub -= c0_;

  // now add linear terms
  for (register int i=0; i<nterms; i++) {

    coeff [i + displacement] = coeff_ [i]; 
    index [i + displacement] = index_ [i];
  }

  // scan arglist for (aux) variables and constants
  for (register int i=0; i<nargs_; i++) {

    expression *curr = arglist_ [i];

    if (curr -> Type () == CONST) {// constant term in sum
      lb -= curr -> Value ();
      ub -= curr -> Value ();
    }
    else {                        // variable
      coeff [++nterms] = 1.; 
      index   [nterms] = curr -> Index ();
    }
  }

  cut -> setRow (nterms + displacement, index, coeff);

  delete [] index;
  delete [] coeff;

  if (lb > -COUENNE_INFINITY) cut -> setLb (lb);
  if (ub <  COUENNE_INFINITY) cut -> setUb (ub);

  cut -> setGloballyValid (); // added only once, it is global
  cs.insert (cut);
  delete cut;
}

/*
 * Name:    linStandardize.cpp
 * Author:  Pietro Belotti
 * Purpose: standardize sum expressions (expr{Sum,Sub,Quad,Group,Opp})
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneProblem.hpp"
#include "exprSum.hpp"
#include "exprSub.hpp"
#include "exprOpp.hpp"
#include "exprMul.hpp"
#include "exprPow.hpp"
#include "exprGroup.hpp"
#include "exprQuad.hpp"

//#define DEBUG

/// analyze sparsity of potential exprQuad/exprGroup and change
/// linear/quadratic maps accordingly, if necessary by adding new
/// auxiliary variables and including them in the linear map
void analyzeSparsity (CouenneProblem *, CouNumber, 
		      std::map <int,                 CouNumber> &,
		      std::map <std::pair <int,int>, CouNumber> &);


/// standardization of linear exprOp's
exprAux *linStandardize (CouenneProblem *p, 
			 bool addAux, 
			 CouNumber c0, 
			 std::map <int,                 CouNumber> &lmap,
 			 std::map <std::pair <int,int>, CouNumber> &qmap) {

  ////////////////////////////////////////////////////////////////////////////////////////

  analyzeSparsity (p, c0, lmap, qmap);

  ////////////////////////////////////////////////////////////////////////////////////////

  int  nq = qmap.size (),     /// data for exprQuad
      *qi = new int [nq+1], 
      *qj = new int [nq+1];
  CouNumber *qc = new CouNumber [nq];

  int  nl = lmap.size(),      /// data for exprGroup
      *li = new int [nl+1];
  CouNumber *lc = new CouNumber [nl];

  // terminate arrays with negative index
  qi [nq] = li [nl] = -1; 

  std::map <int, CouNumber>::iterator lit = lmap.begin (); 

  // fill in arrays for linear part
  for (int i=0; i<nl; i++, lit++) {

    li [i] = lit -> first;
    lc [i] = lit -> second;
  }

  std::map <std::pair <int, int>, CouNumber>::iterator qit = qmap.begin (); 

  // fill in arrays for quadratic part
  for (int i=0; i < nq; i++, qit++) {
    qi [i] = qit -> first. first;
    qj [i] = qit -> first. second;
    qc [i] = qit -> second;
  }

  nl = lmap.size ();
  nq = qmap.size ();

  // particular cases ///////////////////////////////////////////////////////////

  expression *ret;

  // a constant
  if ((nq==0) && (nl==0)) 

    ret = p -> addAuxiliary (new exprConst (c0)); // a constant auxiliary? FIX!

  else if ((nq==0) && (fabs (c0) < COUENNE_EPS) && (nl==1)) { // a linear monomial, cx

    if (fabs (*lc - 1) < COUENNE_EPS) 
      ret    = new exprClone (p -> Var (*li));
    else ret = new exprMul (new exprConst (*lc), new exprClone (p -> Var (*li)));

  } else if ((nl==0) && (fabs (c0) < COUENNE_EPS) && (nq==1)) { 

    // a bilinear/quadratic monomial, cx^2 or cxy

    expression *quad;

    if (*qi == *qj) quad = new exprPow (new exprClone (p -> Var (*qi)), new exprConst (2));
    else            quad = new exprMul (new exprClone (p -> Var (*qi)), 
					new exprClone (p -> Var (*qj)));

    if (fabs (*qc - 1) < COUENNE_EPS) 
      ret    = quad;
    else {
      quad = p -> addAuxiliary (quad);
      ret  = new exprMul (new exprConst (*qc), new exprClone (quad));
    }

  } else {

    // general case ///////////////////////////////////////////////////////////////

    expression **zero = new expression * [1];
    *zero = new exprConst (0.);

    ret = ((nq==0) ? 
      (new exprGroup (c0, li, lc,             zero, 1)) :
      (new exprQuad  (c0, li, lc, qi, qj, qc, zero, 1)));
  }

  delete [] li;
  delete [] lc;
  delete [] qi;
  delete [] qj;
  delete [] qc;

#ifdef DEBUG
  printf ("\nlinstand ==> "); 
  ret -> print (); printf ("\n"); 
  //  ret -> Image () -> print (); printf ("\n");
#endif

  return (addAux ? (p -> addAuxiliary (ret)) : new exprAux (ret));
}

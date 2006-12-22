/*
 * Name:    update_conv.cpp
 * Author:  Pietro Belotti
 * Purpose: update convexification by returning violated convexification cuts 
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>

#include <math.h>

#include <CouenneTypes.h>
#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>

#define SIZE_MALLOC 32768


// reallocates space for buffer buf if n is a multiple of SIZE_MALLOC
// (needs to be called each time a new element is introduced in buf)

extern "C" {
  void Reallocate (void **buf, int n, int size) {

    if (!(n % SIZE_MALLOC) || !(*buf))
      *buf = realloc (*buf, (n + SIZE_MALLOC) * size);
  }
}


// updates problem_'s internal bounds based on user bounds

void CouenneCutGenerator::updateBounds (CouNumber *curlb, 
					CouNumber *curub) {

  // update original variables' bounds

  for (int i = problem_ -> nVars (); i--;) {

    int ind = problem_ -> Var (i) -> Index ();

    problem_ -> Lb (ind) = curlb [ind];
    problem_ -> Ub (ind) = curub [ind];
  }

  // based on the new value for auxiliary variables, update auxiliary bounds

  for (int i = problem_ -> nAuxs (); i--;) {

    exprAux *v = problem_ -> Aux (i);
    int ind = v -> Index ();

    problem_ -> Lb (ind) = (*(v -> Lb ())) ();
    problem_ -> Ub (ind) = (*(v -> Ub ())) ();
  }
}


// Return a set of cuts that tighten convexification.
//
// - x, lb, and ub are the current point, lower-, and upper bound of
// original variables.
//
// The number of generated cuts is returned.

int CouenneCutGenerator::updateConv (CouNumber *curx, 
				     CouNumber *curlb, 
				     CouNumber *curub) {
  // delete all cuts in the pool
  cleanup ();

  // update value of bounds and variables
  expression::update (curx, curlb, curub);

  // At first call, set value of auxiliaries (in subsequent calls,
  // auxiliaries do have their own value)
  if (firstcall_)
    for (int i = problem_ -> nAuxs (); i--;)
      (*(problem_ -> Aux (i))) ();

  updateBounds (curlb, curub);

  /////////////////////////////////// Bound tightening here /////////////////////////////////

  ncuts_ = 0;

  // check for actual violation of coefficients of linear constraints
  // depending on current point
  for (int i = problem_ -> nCons (); i--;) {

    LinearConstraint *con = problem_ -> Con (i);

    CouNumber rhs = (*(con -> Rhs ())) ();

    int nterms = con -> nTerms ();

    CouNumber *coeff = new CouNumber [nterms];

    int *index    = new int [nterms];
    int *conindex = con -> Indices ();
    int  nnz      = 0;

    CouNumber violation = -rhs; 

    CouNumber norm = 0; // l_2 norm of cut

    for (int k=0; k<nterms; k++) {

      CouNumber c = (*(con -> Coeff () [k])) ();

      norm += c*c;

      if (fabs (c) > COUENNE_EPS) {
	coeff [nnz] = c;
	index [nnz] = conindex [k];
	violation += c * curx [conindex [k]];
	nnz++;
      }
    }

    enum con_sign sign = con -> Sign ();

    if      ((sign == COUENNE_GE) && (violation > 0)) violation = 0;
    else if ((sign == COUENNE_LE) && (violation < 0)) violation = 0;

    violation = fabs (violation);

    // checks for violation of exactly one (not both!) of the constraint's bounds
    if (    nnz
	&& (!addviolated_
	    || firstcall_
	    || (violation > COUENNE_EPS))) {

      if (norm > COUENNE_EPS)
	violation /= sqrt (norm);

      OsiRowCut *newcut = new OsiRowCut;

      switch (sign) {
      case COUENNE_EQ: newcut -> setLb (rhs);               newcut -> setUb (rhs);              break;
      case COUENNE_LE: newcut -> setLb (-COUENNE_INFINITY); newcut -> setUb (rhs);              break;
      case COUENNE_GE: newcut -> setLb (rhs);               newcut -> setUb (COUENNE_INFINITY); break;
      default: printf ("Range constraint not used\n"); exit (-1);
      }

      newcut -> setRow (nnz, index, coeff);

      // checks if array pool_ is full and, if so, reallocates space
      // for new entries
      Reallocate ((void **) &pool_, ncuts_, sizeof (OsiRowCut *));

      pool_ [ncuts_++] = newcut;
    }
  }

  pool_ = (OsiRowCut **) realloc (pool_, ncuts_ * sizeof (OsiRowCut *));

  if (firstcall_)
    firstcall_ = false;

  return ncuts_;
}

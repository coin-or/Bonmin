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
/*
extern "C" {
  void Reallocate (void **buf, int n, int size) {

    if (!(n % SIZE_MALLOC) || !(*buf))
      *buf = realloc (*buf, (n + SIZE_MALLOC) * size);
  }
}
*/

// updates problem_'s internal bounds based on user bounds
/*
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
*/

// Return a set of cuts that tighten convexification.
//
// - x, lb, and ub are the current point, lower-, and upper bound of
// original variables.
//
// The number of generated cuts is returned.

int CouenneCutGenerator::updateConv (CouNumber *curx, 
				     CouNumber *curlb, 
				     CouNumber *curub) {

  std::cerr << "--------Update Convexification" << std::endl;

  if (!bonCs_) {

    // This cut generator has been called through updateConv, not
    // through generateCuts, therefore we need to store a
    // OsiSolverInterface somewhere in order to call generateCuts. But
    // first of all, allocate space for bonCs_

    bonCs_ = new OsiCuts;

    // now, create a fake OsiSolverInterface that only needs to
    // contain the value of variables and bounds

    bonOs_ = new OsiClpSolverInterface;

    for (int i=0; i < problem_ -> nVars (); i++)
      bonOs_ -> addCol (0, NULL, NULL, curlb [i], curub [i], 0);
  }
  else {

    // Bonmin is calling this for the second time at least, hence we
    // have to get rid of all cuts contained in the OsiCuts bonCs_
    // before filling it with the new ones.

    for (int i = bonCs_ -> sizeRowCuts (); i--;)
      bonCs_ -> eraseRowCut (i);
  }

  bonOs_ -> setColSolution (curx);

  // ok, now let's just call generateCuts and fill the cuts vector
  // with what we are returned

  generateCuts (*bonOs_, *bonCs_);

  // delete all cuts in the pool
  cleanup ();

  int ncuts = bonCs_ -> sizeRowCuts ();

  if (ncuts) {

    pool_ = (OsiRowCut **) realloc (pool_, ncuts * sizeof (OsiRowCut *));
      for (int i = 0; i<ncuts; i++)
	pool_ [i] = bonCs_ -> rowCutPtr (i);
  }

  printf ("Couenne: %d cuts\n", ncuts);

  return ncuts;

  /*
  // update value of bounds and variables

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
  */
}

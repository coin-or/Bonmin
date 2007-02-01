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


// Return a set of cuts that tighten convexification.
//
// - x, lb, and ub are the current point, lower-, and upper bound of
// original variables.
//
// The number of generated cuts is returned.

int CouenneCutGenerator::updateConv (CouNumber *curx, 
				     CouNumber *curlb, 
				     CouNumber *curub) {
  if (!bonCs_) {

    // This cut generator has been called through updateConv, not
    // through generateCuts, therefore we need to store an
    // OsiSolverInterface and an OsiCuts somewhere in order to call
    // generateCuts. Allocate space for bonCs_

    bonCs_ = new OsiCuts;

    // now, create a fake OsiSolverInterface that only needs to
    // contain the value of variables and bounds

    bonOs_ = new OsiClpSolverInterface;

    int nvars = problem_ -> nVars () + problem_ -> nAuxs ();

    for (int i=0; i < nvars; i++)
      bonOs_ -> addCol (0, NULL, NULL, curlb [i], curub [i], 0);
  }
  else {

    // Bonmin is calling this for the second time at least, hence we
    // have to get rid of all cuts contained in the OsiCuts bonCs_
    // before filling it with the new ones.

    for (int i = bonCs_ -> sizeRowCuts (); i--;)
      bonCs_ -> eraseRowCut (i);

    // update lower and upper bounds
    bonOs_ -> setColLower    (curlb);
    bonOs_ -> setColUpper    (curub);
  }

  bonOs_ -> setColSolution (curx);

  // ok, now let's just call generateCuts and fill the cuts vector
  // with what we are returned

  generateCuts (*bonOs_, *bonCs_);

  // Update pool (used by Bonmin)

  // delete all cuts in the pool
  //  cleanup ();

  ncuts_ = bonCs_ -> sizeRowCuts ();

  // copy cuts into (OsiRowCut ** ) vector for Bonmin's convenience
  if (ncuts_) {

    pool_ = (OsiRowCut **) realloc (pool_, ncuts_ * sizeof (OsiRowCut *));
    for (int i = 0; i<ncuts_; i++)
      pool_ [i] = bonCs_ -> rowCutPtr (i);
  }

  /*

  for (int i=0; i<problem_ -> nVars (); i++) {
    printf ("%4d:  %12.3f [%12.3f %12.3f]\n", 
	    i, problem_ -> X(i), problem_ -> Lb(i), problem_ -> Ub(i));

    //    problem_ -> Var (i) -> print (std::cout);
    //    printf ("\n");
  }

  for (int i=0; i<problem_ -> nAuxs (); i++) {

    int j = i+problem_ -> nVars ();
    printf ("%4d:  %12.3f [%12.3f %12.3f] ", 
	    j, problem_ -> X (j), problem_ -> Lb (j), problem_ -> Ub (j));

    problem_ -> Aux (i) -> print (std::cout);  printf (" = ");
    problem_ -> Aux (i) -> Image () -> print (std::cout);

    expression *lb, *ub;

    problem_ -> Aux (i) -> getBounds (lb, ub);

    printf ("\n L ");
    lb -> print (std::cout);
    printf ("\n U ");
    ub -> print (std::cout);
    printf ("\n");
  }
  */
  return ncuts_;
}

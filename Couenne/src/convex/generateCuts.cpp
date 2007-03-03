/*
 * Name:    generateCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: the generateCuts() method of the convexification class CouenneCutGenerator
 *
 * (C) Pietro Belotti, all rights reserved. 
 * This file is licensed under the Common Public License.
 */

#include <CglCutGenerator.hpp>
#include <CouenneCutGenerator.h>
#include <CouenneProblem.h>


/// a convexifier cut generator

void CouenneCutGenerator::generateCuts (const OsiSolverInterface &si,
					OsiCuts &cs, 
					const CglTreeInfo info) const {
  if (firstcall_) {

    // initialize auxiliary variables and bounds according to originals
    problem_ -> initAuxs (const_cast <CouNumber *> (nlp_ -> getColSolution ()), 
			  const_cast <CouNumber *> (nlp_ -> getColLower    ()),
			  const_cast <CouNumber *> (nlp_ -> getColUpper    ()));

    // OsiSolverInterface is empty yet, no information can be obtained
    // on variables or bounds -- and none is needed since our
    // constructor populated *problem_ with variables and bounds. We
    // only need to update the auxiliary variable and bounds with
    // their current value.

    // add auxiliary variables, unbounded for now
    for (register int i=problem_ -> nAuxs (); i--;)
      const_cast <OsiSolverInterface *> (&si) -> addCol 
	(0, NULL, NULL, - COUENNE_INFINITY, COUENNE_INFINITY, 0);

    // For each auxiliary variable replacing the original (nonlinear)
    // constraints, check if corresponding bounds are violated, and
    // add cut to cs

    int nnlc = problem_ -> nNLCons ();

    for (int i=0; i<nnlc; i++) {

      CouenneConstraint *con = problem_ -> NLCon (i);

      // if there exists violation, add constraint

      OsiRowCut *orc = createCut (0., 0, con -> Body () -> Index (), 1.);

      if (orc) {

	orc -> setLb (con -> Lb () -> Value ());
	orc -> setUb (con -> Ub () -> Value ());
	cs.insert (orc);
	delete orc;
      }
    }
  }

  // first of all, try to tighten the current relaxation by tightening
  // the variables' bounds

  tightenBounds (si);

  // For each auxiliary variable, create cut (or set of cuts) violated
  // by current point and add it to cs

  for (int i=0; i<problem_ -> nAuxs (); i++)
    problem_ -> Aux (i) -> generateCuts (si, cs, this);

  //  if (cs.sizeRowCuts ())
  //    printf ("Couenne: %d convexifier cuts\n", cs.sizeRowCuts ());

  if (firstcall_) {
    firstcall_  = false;
    ntotalcuts_ = nrootcuts_ = cs.sizeRowCuts ();
  }
  else ntotalcuts_ += cs.sizeRowCuts ();
}

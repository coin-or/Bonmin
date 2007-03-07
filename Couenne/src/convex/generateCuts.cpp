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

  //  printf ("/-------------------- GENERATECUTS\n");

  // Contains variables whose bounds have changed due to branching,
  // reduced cost fixing, or bound tightening below. To be used with
  // malloc/realloc/free

  int ncols = problem_ -> nVars () + 
              problem_ -> nAuxs ();

  char *chg_bds = NULL;

  int nchanged = 0; // number of bounds changed;

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
    for (register int i = problem_ -> nVars (), 
	              j = problem_ -> nAuxs (); j--; i++)
      const_cast <OsiSolverInterface *> (&si) -> addCol 
	(0, NULL, NULL, problem_ -> Lb (i), problem_ -> Ub (i), 0);

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
  } else { // equivalent to info.depth > 0

    if (info.pass == 0) // this is the first call in this b&b node
      problem_ -> update (const_cast <CouNumber *> (si. getColSolution ()), 
			  const_cast <CouNumber *> (si. getColLower    ()),
			  const_cast <CouNumber *> (si. getColUpper    ()));

    chg_bds = new char [ncols];

    // fill it with zeros
    for (register int i = ncols; i--;) *chg_bds++ = 0;
    chg_bds -= ncols;

    // not the first call to this procedure, meaning we can be
    // anywhere in the B&B tree. Check, through the auxiliary
    // information, which bounds have changed from the parent node.

    if (info.inTree) {

      OsiBabSolver *auxinfo = dynamic_cast <OsiBabSolver *> (si.getAuxiliaryInfo ());

      if (auxinfo &&
	  (auxinfo -> extraCharacteristics () & 2)) {

	// get previous bounds
	const double * beforeLower = auxinfo -> beforeLower ();
	const double * beforeUpper = auxinfo -> beforeUpper ();

	if (beforeLower && beforeUpper) {

	  // get currentbounds
	  const double * nowLower = si.getColLower();
	  const double * nowUpper = si.getColUpper();

	  for (register int i=0; i < ncols; i++)

	    if ((   nowLower [i] >= beforeLower [i] + COUENNE_EPS)
		|| (nowUpper [i] <= beforeUpper [i] - COUENNE_EPS)) {
	      //	      printf ("x%d bound changed: [%12.4f,%12.4f] -> [%12.4f,%12.4f] ", i,
	      //		      beforeLower [i], beforeUpper [i], 
	      //		      nowLower    [i], nowUpper    [i]);

	      //	      if (i >= problem_ -> nVars ())
	      //		problem_ -> Aux (i - problem_ -> nVars ()) 
	      //		         -> Image () -> print (std::cout);
	      //	      printf ("\n");

	      chg_bds [i] = 1;
	      nchanged++;
	    }

	  //	if (nchanged)
	  //	  printf("%d bounds have changed\n", nchanged);

	} else printf ("WARNING: could not access parent's bounds\n");
      }
    }
  }

  //  for (register int i=0; i < ncols; i++)
  //    printf ("x%3d: [%12.4f,%12.4f]\n", i,
  //	    expression::Lbound (i), 
  //	    expression::Ubound (i));

  // tighten the current relaxation by tightening the variables'
  // bounds

  nchanged = problem_ -> tightenBounds (si, chg_bds);

  //  printf ("::::::::::::::::::::::::::::::::::::::::::::::\n");
  //  for (register int i=0; i < ncols; i++)
  //    printf ("x%3d: [%12.4f,%12.4f]\n", i,
  //	    expression::Lbound (i), 
  //	    expression::Ubound (i));

  //  if (nchanged)
  //    printf("%d bounds tightened\n", nchanged);

  // For each auxiliary variable, create cut (or set of cuts) violated
  // by current point and add it to cs

  // FIXME! Auxiliary has to be re-generated if variables it depends
  // on are in chg_bds, not itself

  //  if (!chg_bds) // this is the first call to generateCuts, we have to
		// generate cuts for all auxiliary variable
    for (int i=0; i < problem_ -> nAuxs (); i++)
      problem_ -> Aux (i) -> generateCuts (si, cs, this);
    /*
  else          // chg_bds contains the indices of the variables whose
		// bounds have changes (a -1 follows the last element)
    for (int i=0, j = problem_ -> nVars (); j < ncols; i++)
      if ((chg_bds [j++]) && 
	  (problem_ -> Aux (i) -> Image () -> Linearity () > LINEAR))
	problem_ -> Aux (i) -> generateCuts (si, cs, this);
    */
  // change tightened bounds through OsiCuts

  if (nchanged) {

    int *indLow = new int [ncols], 
        *indUpp = new int [ncols],
         nLow, nUpp = nLow = 0;

    CouNumber *bndLow = new CouNumber [nchanged],
              *bndUpp = new CouNumber [nchanged];

    const CouNumber 
      *oldLow = si.getColLower (),
      *oldUpp = si.getColUpper (),
      *newLow = problem_ -> Lb (),
      *newUpp = problem_ -> Ub ();

    for (register int i=0; i<ncols; i++)

      if (chg_bds [i]) {

	CouNumber bd;

	if ((bd = newLow [i]) > oldLow [i] + COUENNE_EPS) {
	  indLow [nLow]   = i;
	  bndLow [nLow++] = bd;
	}

	if ((bd = newUpp [i]) < oldUpp [i] - COUENNE_EPS) {
	  indUpp [nUpp]   = i;
	  bndUpp [nUpp++] = bd;
	}
      }

    OsiColCut *cut = new OsiColCut;

    if (cut) {
      cut -> setLbs (nLow, indLow, bndLow);
      cut -> setUbs (nUpp, indUpp, bndUpp);
      cs.insert (cut);
    }

    delete cut;

    delete [] bndLow; delete indLow;
    delete [] bndUpp; delete indUpp;
  }

  if (firstcall_) {
    firstcall_  = false;
    ntotalcuts_ = nrootcuts_ = cs.sizeRowCuts ();
  }
  else ntotalcuts_ += cs.sizeRowCuts ();

  delete [] chg_bds;

  //  printf ("\\________________ generatecuts (%d)\n", cs.sizeRowCuts ());
}

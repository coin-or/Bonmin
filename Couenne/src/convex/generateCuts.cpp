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
#include "BonAuxInfos.hpp"

// fictitious bound for initial unbounded lp relaxations
#define LARGE_BOUND 9.999e12

/// a convexifier cut generator

void CouenneCutGenerator::generateCuts (const OsiSolverInterface &si,
					OsiCuts &cs, 
					const CglTreeInfo info) const {
  
  Bonmin::BabInfo * babInfo = dynamic_cast <Bonmin::BabInfo *> (si.getAuxiliaryInfo ());

  if (babInfo)
    babInfo -> setFeasibleNode ();

  // lift bound on objective auxiliary to avoid overly strict implied
  // bounds

  double now = CoinCpuTime ();

  if (firstcall_) {
    printf ("Couenne:"); 
    fflush (stdout);
  }

  int ncols = problem_ -> nVars () + problem_ -> nAuxs ();

  // This vector contains variables whose bounds have changed due to
  // branching, reduced cost fixing, or bound tightening below. To be
  // used with malloc/realloc/free

  char *chg_bds = new char [ncols];

  // fill it with zeros
  for (register int i = ncols; i--;) 
    *chg_bds++ = 0;
  chg_bds -= ncols;

  if (firstcall_) {

    //////////////////////// FIRST CONVEXIFICATION //////////////////////////////////////

    // initialize auxiliary variables and bounds according to originals
    problem_ -> initAuxs (const_cast <CouNumber *> (nlp_ -> getColSolution ()), 
			  const_cast <CouNumber *> (nlp_ -> getColLower    ()),
			  const_cast <CouNumber *> (nlp_ -> getColUpper    ()));

    // OsiSolverInterface is empty yet, no information can be obtained
    // on variables or bounds -- and none is needed since our
    // constructor populated *problem_ with variables and bounds. We
    // only need to update the auxiliary variables and bounds with
    // their current value.

    // For each auxiliary variable replacing the original (nonlinear)
    // constraints, check if corresponding bounds are violated, and
    // add cut to cs

    int nnlc = problem_ -> nNLCons ();

    for (int i=0; i<nnlc; i++) {

      CouenneConstraint *con = problem_ -> NLCon (i);

      // if there exists violation, add constraint

      int index = con -> Body () -> Index ();

      if (index >= 0) {

	CouNumber l = con -> Lb () -> Value (),	
	          u = con -> Ub () -> Value ();

	// tighten bounds in Couenne's problem representation
	problem_ -> Lb (index) = mymax (l, problem_ -> Lb (index));
	problem_ -> Ub (index) = mymin (u, problem_ -> Ub (index));
      }
    }
  } else { // equivalent to info.depth > 0 || info.pass > 0

    //////////////////////// GET CHANGED BOUNDS DUE TO BRANCHING ////////////////////////

    if (info.pass == 0) // this is the first call in this b&b node
      problem_ -> update (const_cast <CouNumber *> (si. getColSolution ()), 
			  const_cast <CouNumber *> (si. getColLower    ()),
			  const_cast <CouNumber *> (si. getColUpper    ()));

    // not the first call to this procedure, meaning we are anywhere
    // in the B&B tree but at the root node. Check, through the
    // auxiliary information, which bounds have changed from the
    // parent node.

    if (info.inTree) {

      OsiBabSolver *auxinfo = dynamic_cast <OsiBabSolver *> (si.getAuxiliaryInfo ());

      if (auxinfo && (auxinfo -> extraCharacteristics () & 2)) {

	// get previous bounds
	const double * beforeLower = auxinfo -> beforeLower ();
	const double * beforeUpper = auxinfo -> beforeUpper ();

	if (beforeLower && beforeUpper) {

	  // get currentbounds
	  const double * nowLower = si.getColLower();
	  const double * nowUpper = si.getColUpper();

	  for (register int i=0; i < ncols; i++)

	    if ((   nowLower [i] >= beforeLower [i] + COUENNE_EPS)
		|| (nowUpper [i] <= beforeUpper [i] - COUENNE_EPS))
	      chg_bds [i] = 1;

	} else printf ("WARNING: could not access parent's bounds\n");
      }
    }
  }

  //  OsiBabSolver *auxinfo = dynamic_cast <OsiBabSolver *> (si.getAuxiliaryInfo ());

  int objInd = problem_ -> Obj (0) -> Body () -> Index ();
  if (objInd < 0) objInd = 0; 
  /*
  if (babInfo) {
    
    CouNumber primal = babInfo -> mipBound (),
              dual   = babInfo -> bestObjectiveValue ();

    // Bonmin assumes minimization. Hence, primal (dual) is an UPPER
    // (LOWER) bound.

    if (problem_ -> Ub (objInd) > primal) { // update primal bound (MIP)
      problem_ -> Ub (objInd) = primal;
      if (primal < COUENNE_INFINITY - 1) chg_bds [objInd] = 1;
    }

    if (problem_ -> Lb (objInd) < dual) { // update dual bound
      problem_ -> Lb (objInd) = dual;
      if (dual < - COUENNE_INFINITY + 1) chg_bds [objInd] = 1;
    }
  }
  */

  if (BabPtr_) { // update primal bound with best feasible solution object

    int objInd = problem_ -> Obj (0) -> Body () -> Index ();

    if (objInd >= 0) {

      CouNumber bestObj = babInfo -> babPtr() -> model() -> getObjValue();
      
      // Bonmin assumes minimization. Bonmin::Bab::bestObj () should
      // be considered an UPPER bound.

      if (problem_ -> Ub (objInd) > bestObj) {

	problem_ -> Ub (objInd) = bestObj;
	if (bestObj < COUENNE_INFINITY - 1) 
	  chg_bds [objInd] = 1;
      }
    }
  }


  //////////////////////// BOUND TIGHTENING ///////////////////////////////////

  int   ntightened = 0,
      nbwtightened = 0,
      niter = 0, *changed, nchanged;

  do {

    // propagate bounds to auxiliary variables

    ntightened = problem_ -> tightenBounds (chg_bds);

    // implied bounds. Only do it after other bounds have changed,
    // i.e. after branching or if upper bound found

    if (!firstcall_ || chg_bds [objInd]) { 
      //    if (!firstcall_) {

      nbwtightened = problem_ -> impliedBounds (chg_bds);
      
      for (register int i=0; i < ncols; i++) 
	if (expression::Lbound (i) >= expression::Ubound (i) + COUENNE_EPS) {

	  /*printf ("Couenne: infeasible bounds on w_%d [%.12e,%.12e]\n", 
	    i, expression::Lbound (i), expression::Ubound (i));*/

	  /*OsiColCut *infeascut = new OsiColCut;
	  if (infeascut) {
	    double upper = -1., lower = +1.;
	    infeascut -> setLbs (1, &i, &lower);
	    infeascut -> setUbs (1, &i, &upper);
	    cs.insert (infeascut);
	    delete infeascut;
	    }*/

	  //	  infeasNode () = true; // make this node infeasible

	  if (babInfo)
	    babInfo -> setInfeasibleNode ();
	  else printf ("warning, no babinfo\n");

	  goto end_genCuts;
	}
    }
  } while (ntightened && nbwtightened && (niter++ < 10));
  // Not ((ntightened || nbwtightened) && (niter++ < 10)), as we need
  // to repeat only if both bound tighteners had some result.


  //////////////////////// GENERATE CONVEXIFICATION CUTS //////////////////////////////

  // convert sparse chg_bds in something handier
  changed  = (int *) malloc (ncols * sizeof (int));
  nchanged = 0;

  for (register int i=ncols, j=0; i--; j++)
    if (*chg_bds++) {
      *changed++ = j;
      nchanged++;
    }

  delete [] (chg_bds -= ncols);

  changed = (int *) realloc (changed - nchanged, nchanged * sizeof (int));

  // For each auxiliary variable, create cut (or set of cuts) violated
  // by current point and add it to cs

  if (firstcall_)
    for (int i=0, j = problem_ -> nAuxs (); j--; i++)
      problem_ -> Aux (i) -> generateCuts (si, cs, this);

  else { // chg_bds contains the indices of the variables whose bounds
	 // have changes (a -1 follows the last element)

    for (int i=0, j = problem_ -> nAuxs (); j--; i++) {

      expression * image = problem_ -> Aux (i) -> Image ();
      if (   (image -> Linearity () > LINEAR)          // if linear, no need to cut twice
	  && (image -> dependsOn (changed, nchanged))  // if expression does not depend on 
	  )                                            // changed variables, do not cut
	problem_ -> Aux (i) -> generateCuts (si, cs, this);
    }
  }

  if (firstcall_) {

    // set trivial dual bound to objective function, if there is none

    int ind_obj = problem_ -> Obj (0) -> Body () -> Index ();

    if (ind_obj >= 0)
      if (problem_ -> Obj (0) -> Sense () == MINIMIZE) {
	if (problem_ -> Lb (ind_obj) < - LARGE_BOUND)
	  problem_ -> Lb (ind_obj) = - LARGE_BOUND;
      }
      else
	if (problem_ -> Ub (ind_obj) > LARGE_BOUND)
	  problem_ -> Ub (ind_obj) = LARGE_BOUND;
  }

  // change tightened bounds through OsiCuts

  if (nchanged)
    genColCuts (si, cs, nchanged, changed);

  // clean up

  free (changed);

  if (firstcall_) {
    if (cs.sizeRowCuts ())
      printf (" %d initial cuts", cs.sizeRowCuts ());
    printf ("\n");
  }

 end_genCuts:

  if (firstcall_) {
    firstcall_  = false;
    ntotalcuts_ = nrootcuts_ = cs.sizeRowCuts ();
  }
  else ntotalcuts_ += cs.sizeRowCuts ();
					
  septime_ += CoinCpuTime () - now;
}

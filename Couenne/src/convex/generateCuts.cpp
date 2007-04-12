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

// fictitious bound for initial unbounded lp relaxations
#define LARGE_BOUND 9.999e12

/// a convexifier cut generator

void CouenneCutGenerator::generateCuts (const OsiSolverInterface &si,
					OsiCuts &cs, 
					const CglTreeInfo info) const {

  // lift bound on objective auxiliary to avoid overly strict implied
  // bounds

  //  printf ("pass = %d, level = %d, intree = %d, objval = %.12f\n", 
  //	  info.pass, info.level, info.inTree, si.getObjValue());

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

    //for (int i=0; i<15; i++)
    //  printf (":: %3d [%.3e,%.3e]\n", i, problem_ -> Lb (i), problem_ -> Ub (i));

    // initialize auxiliary variables and bounds according to originals
    problem_ -> initAuxs (const_cast <CouNumber *> (nlp_ -> getColSolution ()), 
			  const_cast <CouNumber *> (nlp_ -> getColLower    ()),
			  const_cast <CouNumber *> (nlp_ -> getColUpper    ()));

    //for (int i=0; i<15; i++)
    //  printf ("// %3d [%.3e,%.3e]\n", i, problem_ -> Lb (i), problem_ -> Ub (i));

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

	/*printf ("con %3d [w_%02d]: [%12.3f %12.3f] && [%12.3f %12.3f] --> [%12.3f %12.3f]\n", 
		i, index, 
		l, u, 
		problem_ -> Lb (index), 
		problem_ -> Ub (index),
		mymax (l, problem_ -> Lb (index)),
		mymin (u, problem_ -> Ub (index)));*/

	// tighten bounds in Couenne's problem representation
	problem_ -> Lb (index) = mymax (l, problem_ -> Lb (index));
	problem_ -> Ub (index) = mymin (u, problem_ -> Ub (index));
      }
    }
  } else { // equivalent to info.depth > 0

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
		|| (nowUpper [i] <= beforeUpper [i] - COUENNE_EPS))
	      chg_bds [i] = 1;

	} else printf ("WARNING: could not access parent's bounds\n");
      }
    }
  }

  // update primal bound with best feasible solution object

  if (BabPtr_) {

    int objInd = problem_ -> Obj (0) -> Body () -> Index ();

    if (objInd >= 0) {

      CouNumber bestObj = BabPtr_ -> bestObj();

      if (problem_ -> Obj (0) -> Sense () == MAXIMIZE) { 
	// maximization, bestObj() is a lower bound
	if (problem_ -> Lb (objInd) < bestObj) {
	  //	  printf ("Lower: %.3f", problem_ -> Lb (objInd));
	  problem_ -> Lb (objInd) = bestObj;
	  if (bestObj > - COUENNE_INFINITY + 1) 
	    chg_bds [objInd] = 1;
	  //	  printf (" =-> %.3f\n", problem_ -> Lb (objInd));
	}
      }
      else
	// minimization, bestObj() is an upper bound
	if (problem_ -> Ub (objInd) > bestObj) {
	  //	  printf ("Upper: %.3f", problem_ -> Ub (objInd));
	  problem_ -> Ub (objInd) = bestObj;
	  if (bestObj < COUENNE_INFINITY - 1) 
	    chg_bds [objInd] = 1;
	  //	  printf (" =-> %.3f\n", problem_ -> Ub (objInd));
	}
    }
  }

  //////////////////////// TIGHTEN BOUNDS ///////////////////////////////////

  // tighten the current relaxation by tightening the variables'
  // bounds

  /*printf ("before:");
  for (int i=0; i<problem_ -> nVars () + problem_ -> nAuxs(); i++) {
    printf ("%3d: %+e %+e ", i, problem_ -> Lb (i), problem_ -> Ub (i));
    if (i >= problem_ -> nVars ()) {
      problem_ -> Aux (i-problem_ -> nVars ()) -> print (std::cout);
      printf (" := ");
      problem_ -> Aux (i-problem_ -> nVars ()) -> Image () -> print (std::cout); printf (" [");
      problem_ -> Aux (i-problem_ -> nVars ()) -> Lb () -> print (std::cout); printf (" , ");
      problem_ -> Aux (i-problem_ -> nVars ()) -> Ub () -> print (std::cout); printf ("]");
    }
    printf ("\n");
    }*/
  //(problem_ -> Lb (i) > -COUENNE_INFINITY+1) ? (problem_ -> Lb (i)) : -1e9, 
  //(problem_ -> Ub (i) <  COUENNE_INFINITY-1) ? (problem_ -> Ub (i)) : 1e9);

  int   ntightened = 0, 
      nbwtightened = 0, 
      niter = 0;

  int *changed;
  int  nchanged;

  do {

    // propagate bounds to auxiliary variables

    ntightened = problem_ -> tightenBounds (chg_bds);

    /*for (register int i=0; i < ncols; i++) 
      if (expression::Lbound (i) >= expression::Ubound (i) + COUENNE_EPS)
	printf ("Couenne: infeasible bounds on w_%d [%.12e,%.12e]\n", 
	i, expression::Lbound (i), expression::Ubound (i));*/

    // implied bounds. Only do it after other bounds have changed,
    // i.e. after branching

    if (!firstcall_) {

      nbwtightened = problem_ -> impliedBounds (chg_bds);
      
      /*if (ntightened || nbwtightened) printf ("%d, %d bounds tightened\n", 
	ntightened, nbwtightened);*/

      for (register int i=0; i < ncols; i++) 
	if (expression::Lbound (i) >= expression::Ubound (i) + COUENNE_EPS) {

	  /*printf ("Couenne: infeasible bounds on w_%d [%.12e,%.12e]\n", 
	    i, expression::Lbound (i), expression::Ubound (i));*/

	  OsiColCut *infeascut = new OsiColCut;

	  if (infeascut) {

	    double upper = -1, lower = +1; // want to make this node infeasible

	    infeascut -> setLbs (1, &i, &lower);
	    infeascut -> setUbs (1, &i, &upper);

	    cs.insert (infeascut);
	    delete infeascut;

	    // avoid bookkeeping of infeasible flag
	    goto end_genCuts;
	  }
	}
    }
  } while (ntightened && nbwtightened && (niter++ < 10));
  // as opposed to ((ntightened || nbwtightened) && (niter++ < 10));


  //////////////////////// GENERATE CONVEXIFICATION CUTS //////////////////////////////

  // first of all, convert sparse chg_bds in something handier

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
      /*
      if (   (image -> Linearity () > LINEAR)          // if linear, no need to cut twice
	  && (image -> dependsOn (changed, nchanged))  // if expression does not depend on 
	  )                                            // changed variables, do not cut
      */
	problem_ -> Aux (i) -> generateCuts (si, cs, this);
    }
  }

  //////////////////////// GENERATE OsiColCut FOR SHRUNKEN BOUNDS /////////////////////

  if (firstcall_) {

    // set trivial dual bound to objective function

    int ind_obj = problem_ -> Obj (0) -> Body () -> Index ();

    if (ind_obj >= 0)
      if (problem_ -> Obj (0) -> Sense () == MINIMIZE)
	problem_ -> Lb (ind_obj) = - LARGE_BOUND;
      else
	problem_ -> Ub (ind_obj) =   LARGE_BOUND;
  }

  // change tightened bounds through OsiCuts

  if (nchanged) {

    // indices for OsiColCut
    int *indLow = new int [ncols], 
        *indUpp = new int [ncols],
         nLow, nUpp = nLow = 0;

    // values fo OsiColCut
    CouNumber *bndLow = new CouNumber [ncols],
              *bndUpp = new CouNumber [ncols];

    const CouNumber 
      *oldLow = si.getColLower (), // old bounds
      *oldUpp = si.getColUpper (),
      *newLow = problem_ -> Lb (), // changed bounds
      *newUpp = problem_ -> Ub ();

    // check all changed bounds
    for (int i=0; i<nchanged; i++) {

      register int index = changed [i];

      CouNumber bd;

      if (firstcall_ || ((bd = newLow [index]) > oldLow [index] + COUENNE_EPS)) { // lower
	indLow [nLow]   = index;
	bndLow [nLow++] = bd;
      }

      if (firstcall_ || ((bd = newUpp [index]) < oldUpp [index] - COUENNE_EPS)) { // upper
	indUpp [nUpp]   = index;
	bndUpp [nUpp++] = bd;
      }
    }

    // ok, now create Column Cut

    if (nUpp || nLow) {

      OsiColCut *cut = new OsiColCut;

      if (cut) {
	cut -> setLbs (nLow, indLow, bndLow);
	cut -> setUbs (nUpp, indUpp, bndUpp);

	cs.insert (cut);
	delete cut;
      }
    }

    delete [] bndLow; delete [] indLow;
    delete [] bndUpp; delete [] indUpp;
  }

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
					
  /*printf ("after: ");
  for (int i=0; i<problem_ -> nVars () + problem_ -> nAuxs(); i++)
    printf ("%.1f %.1f|",
	    (problem_ -> Lb (i) > -COUENNE_INFINITY+1) ? (problem_ -> Lb (i)) : -1e9, 
	    (problem_ -> Ub (i) <  COUENNE_INFINITY-1) ? (problem_ -> Ub (i)) : 1e9);
	    printf ("\n");*/

  septime_ += CoinCpuTime () - now;
}

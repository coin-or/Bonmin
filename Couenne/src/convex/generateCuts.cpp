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


// a convexifier cut generator

void CouenneCutGenerator::generateCuts (const OsiSolverInterface &si, 
					OsiCuts &cs, 
					const CglTreeInfo info) const {

  if (firstcall_) {

    // OsiSolverInterface is empty yet, no information can be obtained
    // on variables or bounds -- and none is needed since our
    // constructor populated *problem_ with variables and bounds. We
    // only need to 

    // For each auxiliary variable replacing the original constraints,
    // check if corresponding bounds are violated, and add cut to cs

    int nnlc = problem_ -> nNLCons ();

    for (int i=0; i<nnlc; i++) {

      CouenneConstraint *con = problem_ -> NLCon (i);

      // for constraint lb <= w <= ub, compute actual values of lb, ub

      CouNumber lb = con -> Lb () -> Value ();
      CouNumber ub = con -> Ub () -> Value ();

      // if there exists violation, add constraint

      OsiRowCut *orc = createCut (0., 0, con -> Body () -> Index (), 1.);

      /*
      OsiRowCut *orc   = new OsiRowCut;
      CouNumber *coeff = new CouNumber [1];
      int       *index = new int       [1];

      coeff [0] = 1;
      index [0] = con -> Body () -> Index ();

      orc -> setRow (1, index, coeff);
      */

      if (orc) {
	if (lb > - COUENNE_INFINITY + 1) orc -> setLb (lb);
	if (ub <   COUENNE_INFINITY - 1) orc -> setUb (ub);
	
	cs.insert (orc);

	delete orc;
      }
      //      delete [] coeff;
      //      delete [] index;
    }
  }
  else {

    // Retrieve, from si, value and bounds of all variables, if not
    // firstcall, otherwise only those of the original ones Update
    // expression structure with x, l, u

    CouNumber 
      *xc = const_cast <CouNumber *> (si.getColSolution ()), 
      *lc = const_cast <CouNumber *> (si.getColLower    ()),
      *uc = const_cast <CouNumber *> (si.getColUpper    ());

    // update now all variables and bounds

    problem_ -> update (xc, lc, uc);

    // update bounding box (which may depend on the original
    // variables' box) for the variables whose bound is looser. Here,
    // newly enforced branching rules may change dependent auxiliary
    // variables' bounds, in a recursive way. Hence we need to repeat
    // the propagation step as long as at least one bound is modified.

    bool found_one;

    do {

      found_one = false;

      int naux = problem_ -> nAuxs ();

      // check all auxiliary variables for changes in their upper,
      // lower bound, depending on the bound changes of the variables
      // they depend on

      for (register int i = problem_ -> nVars (), j=0; 
	   j < naux; j++) {
    
	CouNumber ll = (*(problem_ -> Aux (j) -> Lb ())) ();
	CouNumber uu = (*(problem_ -> Aux (j) -> Ub ())) ();

	// check if lower bound got higher    
	if (ll > lc [i+j]) {
	  const_cast <OsiSolverInterface *> (&si) -> setColLower (i+j, ll);
	  lc [i+j] = ll;
	  found_one = true;
	}

	// check if upper bound got lower
	if (uu < uc [i+j]) {
	  const_cast <OsiSolverInterface *> (&si) -> setColUpper (i+j, uu);
	  uc [i+j] = uu;
	  found_one = true;
	}

	expression::update (xc, lc, uc);
      }

    } while (found_one); // repeat as long as at least one bound changed

    // update again 
    problem_ -> update (xc, lc, uc);
  }

  // For each auxiliary variable, create cut (or set of cuts) violated
  // by current point and add it to cs

  for (int i=0; i<problem_ -> nAuxs (); i++)
    problem_ -> Aux (i) -> generateCuts (si, cs, this);

  // end of generateCuts

  if (cs.sizeRowCuts ())
    printf ("Couenne: %d convexifier cuts\n", cs.sizeRowCuts ());

  if (firstcall_) 
    firstcall_ = false;
}

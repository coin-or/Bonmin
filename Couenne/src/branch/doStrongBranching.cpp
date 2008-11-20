/*
 * Name:    doStrongBranching.cpp
 * Authors: Andreas Waechter, IBM Corp.
 *          Pietro Belotti, CMU
 * Purpose: actual strong branching method
 *
 * (C) Carnegie-Mellon University, 2008.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinTime.hpp"
#include "CouenneChooseStrong.hpp"
#include "CouenneProblem.hpp"
#include "CouenneBranchingObject.hpp"

/// compute Euclidean distance between two points (most likely LP solutions)
/// l_2 norm by default, but can change it by fourth parameter
double distance (const double *p1, const double *p2, int size, double k=2.) {

  double 
    result = 0.,
    element;

  if (k == 2.) // a bit faster, probably

    while (size--) {
      element = *p1++ - *p2++;
      result += element * element;
    }

  else

    while (size--) {
      element = *p1++ - *p2++;
      result += pow (element, k);
    }

  return pow (result, 1. / k);
}


namespace Bonmin {

  /**  This is a utility function which does strong branching on
       a list of objects and stores the results in OsiHotInfo.objects.
       On entry the object sequence is stored in the OsiHotInfo object
       and maybe more.
       It returns -

      -1 - one branch was infeasible both ways
       0 - all inspected - nothing can be fixed
       1 - all inspected - some can be fixed (returnCriterion==0)
       2 - may be returning early - one can be fixed (last one done) (returnCriterion==1) 
       3 - returning because max time
  */
  int CouenneChooseStrong::doStrongBranching (OsiSolverInterface * solver, 
					      OsiBranchingInformation *info,
					      int numberToDo, int returnCriterion)
  {

    jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, 
		      "\n-\n------- CCS: trying %d objects:\n", numberToDo);

    int numberColumns = solver -> getNumCols ();

    solver -> markHotStart (); // save current LP point

    const double
      *lower = info -> lower_,
      *upper = info -> upper_;

    double 
      *saveLower = CoinCopyOfArray (info -> lower_, numberColumns),
      *saveUpper = CoinCopyOfArray (info -> upper_, numberColumns),
      *Lower0    = CoinCopyOfArray (info -> lower_, numberColumns), // delete afterwards
      *Upper0    = CoinCopyOfArray (info -> upper_, numberColumns),
      *oldLower  = new double [numberColumns],
      *oldUpper  = new double [numberColumns],
      *lpSol     = NULL, 
       timeStart = CoinCpuTime ();

    // LP solution for distance
    if (pseudoUpdateLP_) 
      lpSol = CoinCopyOfArray (info -> solution_, numberColumns);

    // provide Couenne problem with point/bounds contained in info
    problem_ -> domain () -> push
      (problem_ -> nVars (),
       info -> solution_,
       info -> lower_,
       info -> upper_);

    int returnCode = 0, iDo = 0;

    for (;iDo < numberToDo; iDo++) {

      HotInfo * result = results_ () + iDo; // retrieve i-th object to test

      CouenneObject *CouObj = dynamic_cast <CouenneObject *>
	(solver_ -> objects () [result -> whichObject ()]);

      // For now just 2 way
      OsiBranchingObject * branch = result -> branchingObject ();
      assert (branch->numberBranches()==2);

      CouenneBranchingObject *cb = dynamic_cast <CouenneBranchingObject *> (branch);
      if (cb) cb -> setSimulate (true);

      /* Try the first direction.  Each subsequent call to branch()
	 performs the specified branch and advances the branch object
	 state to the next branch alternative.) */

      int 
	status0 = -1, 
	status1 = -1;

      OsiSolverInterface * thisSolver = solver; 

      // DOWN DIRECTION ///////////////////////////////////////////////////////

      if (branch -> boundBranch ()) { // a (variable) bound branch

        if (branch -> branch (solver) > COUENNE_INFINITY) // branch is infeasible
	  result -> setDownStatus (status0 = 1);

	else { // branch is feasible, solve and compare

	  solver -> solveFromHotStart ();
	  if (pseudoUpdateLP_ && CouObj && solver -> isProvenOptimal ()) {
	    CouNumber dist = distance (lpSol, solver -> getColSolution (), numberColumns);
	    if (dist > COUENNE_EPS)
	      CouObj -> setEstimate (dist, 0);
	  }
	}

      } else {                       // some more complex branch, have to clone solver

        // adding cuts or something 
        thisSolver = solver -> clone ();

        if (branch -> branch (thisSolver) > COUENNE_INFINITY)
	  result -> setDownStatus (status0 = 1);

	else { // set hot start iterations
	  int limit;
	  thisSolver -> getIntParam (OsiMaxNumIterationHotStart, limit);
	  thisSolver -> setIntParam (OsiMaxNumIteration,         limit); 

	  thisSolver -> resolve ();
	  if (pseudoUpdateLP_ && CouObj && thisSolver -> isProvenOptimal ()) {
	    CouNumber dist = distance (lpSol, thisSolver -> getColSolution (), numberColumns);
	    if (dist > COUENNE_EPS)
	      CouObj -> setEstimate (dist, 0);
	    //CouObj -> setEstimate (distance (lpSol, thisSolver->getColSolution (),numberColumns), 0);
	  }
	}
      }

      // can check if we got solution
      // status is 0 finished, 1 infeasible and 2 unfinished and 3 is solution

      // only update information if this branch is feasible
      if (status0 < 0)
	status0 = result -> updateInformation (thisSolver, info, this);

      numberStrongIterations_ += thisSolver -> getIterationCount();

      if ((status0 == 3) && (trustStrongForSolution_)) {
        // new solution already saved
  	info -> cutoff_ = goodObjectiveValue_;
	problem_ -> setCutOff (goodObjectiveValue_);
  	status0 = 0;
      }

      if (solver != thisSolver)
        delete thisSolver;

      // save current bounds as tightened by the down branch; will be
      // used below to update global bounding box in solver
      CoinCopyN (problem_ -> Lb (), numberColumns, oldLower);
      CoinCopyN (problem_ -> Ub (), numberColumns, oldUpper);

      // Restore pre-left-branch bounds in solver
      for (int j=0; j<numberColumns; j++) {

        if (saveLower [j] != lower [j]) solver -> setColLower (j, saveLower [j]);
        if (saveUpper [j] != upper [j]) solver -> setColUpper (j, saveUpper [j]);
      }

      // UP DIRECTION ///////////////////////////////////////////////////////

      thisSolver = solver; 

      if (branch -> boundBranch ()) { // (variable) bound branch 

        if (branch -> branch (solver) > COUENNE_INFINITY)
	  result -> setUpStatus (status1 = 1);

        else {
	  solver -> solveFromHotStart ();
	  if (pseudoUpdateLP_ && CouObj && solver -> isProvenOptimal ()) {
	    CouNumber dist = distance (lpSol, solver -> getColSolution (), numberColumns);
	    if (dist > COUENNE_EPS)
	      CouObj -> setEstimate (dist, 1);
	    //CouObj -> setEstimate (distance (lpSol, solver -> getColSolution (), numberColumns), 1);
	  }
	}
      } else {                     // some more complex branch, have to clone solver
        // adding cuts or something 
        thisSolver = solver -> clone ();

        if (branch -> branch (thisSolver) > COUENNE_INFINITY)
	  result -> setUpStatus (status1 = 1);

	else {
	  // set hot start iterations
	  int limit;
	  thisSolver -> getIntParam (OsiMaxNumIterationHotStart, limit);
	  thisSolver -> setIntParam (OsiMaxNumIteration,         limit); 

	  thisSolver -> resolve();
	  if (pseudoUpdateLP_ && CouObj && thisSolver -> isProvenOptimal ()) {
	    CouNumber dist = distance (lpSol, thisSolver -> getColSolution (), numberColumns);
	    if (dist > COUENNE_EPS)
	      CouObj -> setEstimate (dist, 1);
	    //CouObj -> setEstimate (distance (lpSol, thisSolver->getColSolution (),numberColumns), 1);
	  }
	}
      }

      // can check if we got solution
      // status is 0 finished, 1 infeasible and 2 unfinished and 3 is solution

      // only update information if this branch is feasible
      if (status1 < 0)
	status1 = result -> updateInformation (thisSolver, info, this);

      numberStrongDone_++;
      numberStrongIterations_ += thisSolver->getIterationCount();

      if ((status1 == 3) && (trustStrongForSolution_)) {
        // new solution already saved
	info -> cutoff_ = goodObjectiveValue_;
	problem_ -> setCutOff (goodObjectiveValue_);
	status1 = 0;
      }

      jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, "-------\n");

      if (cb) cb -> setSimulate (false);

      /////////////////////////////////////////////////////////////////////////////

      if (solver != thisSolver)
        delete thisSolver;

      bool tightened = false;

      t_chg_bounds *chg_bds = new t_chg_bounds [numberColumns];

      // extend problem_'s bounding box to include downbranch's tightened
      for (int j=0; j<numberColumns; j++) {

        if (oldLower [j] < problem_ -> Lb (j)) problem_ -> Lb (j) = oldLower [j];
        if (oldUpper [j] > problem_ -> Ub (j)) problem_ -> Ub (j) = oldUpper [j];

	if (problem_ -> Lb (j) > lower [j] + COUENNE_EPS) {
	  chg_bds [j].setLower (t_chg_bounds::CHANGED);
	  tightened = true;
	}

	if (problem_ -> Ub (j) < upper [j] - COUENNE_EPS) {
	  chg_bds [j].setUpper (t_chg_bounds::CHANGED);
	  tightened = true;
	}
      }

      if (tightened &&                     // have tighter bounds
	  (problem_ -> doFBBT ()) &&       // selected FBBT
	  !(problem_ -> btCore (chg_bds))) // tighten again on root

	status0 = status1 = 1;	           // if returns false, problem is infeasible

      delete [] chg_bds;

      // create union of bounding box from both branching directions
      for (int j=0; j<numberColumns; j++) {

        if (oldLower [j] < problem_ -> Lb (j)) problem_ -> Lb (j) = oldLower [j];
        if (oldUpper [j] > problem_ -> Ub (j)) problem_ -> Ub (j) = oldUpper [j];
      }

      // set new bounding box as the possibly tightened one (a subset
      // of the initial)
      for (int j=0; j<numberColumns; j++) {

        solver -> setColLower (j, saveLower [j] = problem_ -> Lb (j));
        solver -> setColUpper (j, saveUpper [j] = problem_ -> Ub (j));
      }

      /*
        End of evaluation for this candidate object. Possibilities are:

        * Both sides below cutoff; this variable is a candidate for
          branching.

        * Both sides infeasible or above the objective cutoff: no
          further action here. Break from the evaluation loop and
          assume the node will be purged by the caller.

        * One side feasible and below cutoff: Install the branch
          (i.e., fix the variable). Possibly break from the evaluation
          loop and assume the node will be reoptimised by the caller.
      */

      if (status0 == 1 && 
	  status1 == 1) { // infeasible
        returnCode=-1;
        break; // exit loop
      } else if (status0==1 || status1==1) {
        numberStrongFixed_++;
        if (!returnCriterion) {
	  returnCode=1;
        } else {
	  returnCode=2;
	  break;
        }
      }

      bool hitMaxTime = ( CoinCpuTime()-timeStart > info->timeRemaining_);
      if (hitMaxTime) {
        returnCode=3;
        break;
      }
    } // end loop /***********************************/
  

    if (jnlst_ -> ProduceOutput (J_DETAILED, J_BRANCHING)) {
      printf ("tightened bounds: ");
      // create union of bounding box from both branching directions
      for (int j=0; j<numberColumns; j++) {
      
	if (problem_ -> Lb (j) > Lower0 [j]) printf ("l%d (%g-->%g) ", j,Lower0[j], problem_->Lb (j));
	if (problem_ -> Ub (j) < Upper0 [j]) printf ("u%d (%g-->%g) ", j,Upper0[j], problem_->Ub (j));
      }
    }

    delete [] Lower0;
    delete [] Upper0;

    problem_ -> domain () -> pop (); // discard current point/bounds from problem

    delete [] lpSol;

    jnlst_ -> Printf (J_ITERSUMMARY, J_BRANCHING, "----------------------done\n\n\n");

    if (iDo < numberToDo) iDo++; // exited due to infeasibility
    assert (iDo <= (int) results_.size());
    results_.resize (iDo);

    delete [] oldLower;
    delete [] oldUpper;

    delete [] saveLower;
    delete [] saveUpper;

    solver -> unmarkHotStart ();     // Delete the snapshot

    return returnCode;
  }
}

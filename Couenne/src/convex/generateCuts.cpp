/*
 * Name:    generateCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: the generateCuts() method of the convexification class CouenneCutGenerator
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "BonAuxInfos.hpp"
#include "CglCutGenerator.hpp"

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "CouenneSolverInterface.hpp"

// exception
#define INFEASIBLE 1

// fictitious bound for initial unbounded lp relaxations
#define LARGE_BOUND 9.999e12

#define LARGE_TOL (LARGE_BOUND / 1e6)

// set and lift bound for auxiliary variable associated with objective
// function
void fictitiousBound (OsiCuts &cs,
		      CouenneProblem *p, 
		      bool action) {     // true before convexifying, false afterwards

  // set trivial dual bound to objective function, if there is none

  int ind_obj = p -> Obj (0) -> Body () -> Index ();

  if (ind_obj < 0) return;

  // we have a single variable objective function

  int sense = (p -> Obj (0) -> Sense () == MINIMIZE) ? -1 : 1;

  if (action)
    if (sense<0) {if (p -> Lb (ind_obj) < - LARGE_BOUND) p -> Lb (ind_obj) = - LARGE_BOUND;}
    else         {if (p -> Ub (ind_obj) >   LARGE_BOUND) p -> Ub (ind_obj) =   LARGE_BOUND;}
  else
    if (sense>0) {if (fabs (p->Ub(ind_obj)-LARGE_BOUND)<LARGE_TOL) p->Ub(ind_obj) = COUENNE_INFINITY;}
    else         {if (fabs (p->Lb(ind_obj)+LARGE_BOUND)<LARGE_TOL) p->Lb(ind_obj) =-COUENNE_INFINITY;}
}


// translate changed bound sparse array into a dense one
void sparse2dense (int ncols, t_chg_bounds *chg_bds, int *&changed, int &nchanged) {

  // convert sparse chg_bds in something handier
  // AW: replacd "malloc" here by "realloc"; otherwise this is a memory leak
  //     In general, I don't think it is worth to do a realloc here,
  //     it is probably more expensive than not using it.  The memory
  //     is free anyway when generateCuts is left
  changed  = (int *) realloc (changed, ncols * sizeof (int));
  nchanged = 0;

  for (register int i=ncols, j=0; i--; j++, chg_bds++)
    if (chg_bds -> lower() != t_chg_bounds::UNCHANGED ||
	chg_bds -> upper() != t_chg_bounds::UNCHANGED ) {
      *changed++ = j;
      nchanged++;
    }

  changed -= nchanged;
  //changed = (int *) realloc (changed, nchanged * sizeof (int));
}


/// get new bounds from parents' bounds + branching rules
void updateBranchInfo (const OsiSolverInterface &si, CouenneProblem *p, 
		       t_chg_bounds *chg, const CglTreeInfo &info);

/// a convexifier cut generator

void CouenneCutGenerator::generateCuts (const OsiSolverInterface &si,
					OsiCuts &cs, 
					const CglTreeInfo info) const {
  int nInitCuts = cs.sizeRowCuts ();

  /*static int count = 0;
  char fname [20];
  sprintf (fname, "relax_%d", count++);
  si.writeLp (fname);
  printf ("writing %s\n", fname);*/

  jnlst_ -> Printf (J_DETAILED, J_CONVEXIFYING,
		    "generateCuts: level = %d, pass = %d, intree = %d\n",
		    info.level, info.pass, info.inTree);

  Bonmin::BabInfo * babInfo = dynamic_cast <Bonmin::BabInfo *> (si.getAuxiliaryInfo ());

  if (babInfo)
    babInfo -> setFeasibleNode ();

  double now   = CoinCpuTime ();
  int    ncols = problem_ -> nVars ();

  problem_ -> installCutOff ();

  // This vector contains variables whose bounds have changed due to
  // branching, reduced cost fixing, or bound tightening below. To be
  // used with malloc/realloc/free

  t_chg_bounds *chg_bds = new t_chg_bounds [ncols];

  if (firstcall_) {

    // First convexification //////////////////////////////////////

    // OsiSolverInterface is empty yet, no information can be obtained
    // on variables or bounds -- and none is needed since our
    // constructor populated *problem_ with variables and bounds. We
    // only need to update the auxiliary variables and bounds with
    // their current value.

    // For each auxiliary variable replacing the original (nonlinear)
    // constraints, check if corresponding bounds are violated, and
    // add cut to cs

    int nnlc = problem_ -> nCons ();

    for (int i=0; i<nnlc; i++) {

      // for each constraint
      CouenneConstraint *con = problem_ -> Con (i);

      // (which has an aux as its body)
      int index = con -> Body () -> Index ();

      if ((index >= 0) && (con -> Body () -> Type () == AUX)) {

	// get the auxiliary that is at the lhs
	exprAux *conaux = dynamic_cast <exprAux *> (problem_ -> Var (index));

	if (conaux &&
	    (conaux -> Image ()) && 
	    (conaux -> Image () -> Linearity () <= LINEAR)) {

	  // the auxiliary w of constraint w <= b is associated with a
	  // linear expression w = ax: add constraint ax <= b
	  conaux -> Image () -> generateCuts (conaux, si, cs, this, chg_bds, 
					      conaux -> Index (), 
					      (*(con -> Lb ())) (), 
					      (*(con -> Ub ())) ());

	  // take it from the list of the variables to be linearized
	  // 
	  // DO NOT decrease multiplicity. Even if it is a linear
	  // term, its bounds can still be used in implied bounds
	  //
	  //conaux -> decreaseMult (); // !!!
	}

	// also, add constraint w <= b

	// if there exists violation, add constraint
	CouNumber l = con -> Lb () -> Value (),	
	          u = con -> Ub () -> Value ();

	//printf ("constraint %d: [%g,%g]", i, l, u); con -> print ();

	// tighten bounds in Couenne's problem representation
	problem_ -> Lb (index) = CoinMax (l, problem_ -> Lb (index));
	problem_ -> Ub (index) = CoinMin (u, problem_ -> Ub (index));

      } else { // body is more than just a variable, but it should be
	       // linear. If so, generate equivalent linear cut

	assert (false);	// TODO
      }
    }

    if (jnlst_ -> ProduceOutput (J_DETAILED, J_CONVEXIFYING)) {
      if (cs.sizeRowCuts ()) {
	jnlst_ -> Printf (J_DETAILED, J_CONVEXIFYING,"Couenne: constraint row cuts\n");
	for (int i=0; i<cs.sizeRowCuts (); i++) 
	  cs.rowCutPtr (i) -> print ();
      }
      if (cs.sizeColCuts ()) {
	jnlst_ -> Printf (J_DETAILED, J_CONVEXIFYING,"Couenne: constraint col cuts\n");
	for (int i=0; i<cs.sizeColCuts (); i++) 
	  cs.colCutPtr (i) -> print ();
      }
    }
  } else {

    // use new optimum as lower bound for variable associated w/objective
    int indobj = problem_ -> Obj (0) -> Body () -> Index ();

    //assert (indobj >= 0);

    /*CouNumber save_obj_primal = 
      (problem_ -> Obj (0) -> Sense () == MINIMIZE) ? 
      problem_ -> Ub (indobj) : 
      problem_ -> Lb (indobj);*/

    // transmit solution from OsiSolverInterface to problem
    problem_ -> domain () -> push 
      (problem_ -> nVars (),
       si. getColSolution (), 
       si. getColLower    (),
       si. getColUpper    ());

    /*    ((problem_ -> Obj (0) -> Sense () == MINIMIZE) ? 
     problem_ -> Ub (indobj) : 
     problem_ -> Lb (indobj)) = save_obj_primal;*/

    if (indobj >= 0) {

      /*if (problem_ -> Obj (0) -> Sense () == MINIMIZE) {
	const_cast <OsiSolverInterface *> (&si) ->
	  setColUpper (indobj, problem_ -> getCutOff ());
	problem_ -> domain () -> ub (indobj) = problem_ -> getCutOff ();
      } else {
	const_cast <OsiSolverInterface *> (&si) ->
	  setColLower (indobj, problem_ -> getCutOff ());
	problem_ -> domain () -> lb (indobj) = problem_ -> getCutOff ();
	}*/

      // Use current value of objvalue's x as a lower bound for bound
      // tightening
      double lp_bound = problem_ -> domain () -> x (indobj);

      if (problem_ -> Obj (0) -> Sense () == MINIMIZE) 
	   {if (lp_bound > problem_ -> Lb (indobj)) problem_ -> Lb (indobj) = lp_bound;}
      else {if (lp_bound < problem_ -> Ub (indobj)) problem_ -> Ub (indobj) = lp_bound;}
    }

    updateBranchInfo (si, problem_, chg_bds, info); // info.depth >= 0 || info.pass >= 0
  }

  fictitiousBound (cs, problem_, false);

  problem_ -> installCutOff ();

  if (jnlst_ -> ProduceOutput (J_VECTOR, J_CONVEXIFYING)) {
    jnlst_ -> Printf(J_VECTOR, J_CONVEXIFYING,"== generateCuts: point to cut =============\n");
    for (int i = 0; i < problem_ -> nVars (); i++)
      if (problem_ -> Var (i) -> Multiplicity () > 0)
	jnlst_->Printf(J_VECTOR, J_CONVEXIFYING,
		       "%4d %+20.8g [%+20.8g,%+20.8g] --- %+20.8g [%+20.8g,%+20.8g] (%+20.8g)\n", i,
		       problem_ -> X  (i),
		       problem_ -> Lb (i),
		       problem_ -> Ub (i),
		       si.getColSolution  () [i],
		       si.getColLower  () [i],
		       si.getColUpper  () [i],
		       problem_ -> bestSol () ? problem_ -> bestSol () [i] : 0.);
    jnlst_->Printf(J_VECTOR, J_CONVEXIFYING,"=============================\n");
  }

  int *changed = NULL, nchanged;

  // Bound tightening ///////////////////////////////////////////

  // do bound tightening only at first pass of cutting plane in a node
  // of BB tree (info.pass == 0) or if first call (creation of RLT,
  // info.pass == -1)

  problem_ -> installCutOff ();

  try {

    // Bound tightening ////////////////////////////////////

    // FBBT
    if (problem_ -> doFBBT () && (info.pass <= 0) &&
	(! (problem_ -> boundTightening (chg_bds, babInfo))))
      throw INFEASIBLE;

    // OBBT
    if (!firstcall_ && // no obbt if first call (there is no LP to work with)
	problem_ -> obbt (this, si, cs, info, babInfo, nchanged, changed, chg_bds) < 0)
      throw INFEASIBLE;

    // Reduced Cost BT
    if (problem_ -> doFBBT () && !firstcall_)
      problem_ -> redCostBT (&si, chg_bds, babInfo);

    if ((problem_ -> doFBBT () ||
	 problem_ -> doOBBT () ||
	 problem_ -> doABT  ()) &&
	(jnlst_ -> ProduceOutput (J_VECTOR, J_CONVEXIFYING))) {
      jnlst_ -> Printf(J_VECTOR, J_CONVEXIFYING,"== after bt =============\n");
      for (int i = 0; i < problem_ -> nVars (); i++)
	if (problem_ -> Var (i) -> Multiplicity () > 0)
	  jnlst_->Printf(J_VECTOR, J_CONVEXIFYING,"%4d %+20.8g [%+20.8g,%+20.8g]\n", i,
			 problem_ -> X  (i),
			 problem_ -> Lb (i),
			 problem_ -> Ub (i));
      jnlst_->Printf(J_VECTOR, J_CONVEXIFYING,"=============================\n");
    }

    // Generate convexification cuts //////////////////////////////

    sparse2dense (ncols, chg_bds, changed, nchanged);

    double *nlpSol;

    //--------------------------------------------

    if (babInfo && ((nlpSol = const_cast <double *> (babInfo -> nlpSolution ())))) {

      // Aggressive Bound Tightening ////////////////////////////////

      int logAbtLev = problem_ -> logAbtLev ();

      if (problem_ -> doABT () &&           // flag is checked, AND
	  ((logAbtLev != 0) ||                // (parameter is nonzero OR
	   (info.level == 0)) &&              //  we are at root node), AND
	  (info.pass == 0) &&               // at first round of cuts, AND 
	  ((logAbtLev < 0) ||                 // (logAbtLev = -1, OR
	   (info.level <= logAbtLev) ||       //  depth is lower than COU_OBBT_CUTOFF_LEVEL, OR
	   (CoinDrand48 () <                  //  probability inversely proportional to the level)
	    pow (2., (double) logAbtLev - (info.level + 1))))) {

	jnlst_ -> Printf(J_DETAILED, J_CONVEXIFYING,"  performing ABT\n");
	if (! (problem_ -> aggressiveBT (nlp_, chg_bds, babInfo)))
	  throw INFEASIBLE;

	sparse2dense (ncols, chg_bds, changed, nchanged);
      }

      // obtain solution just found by nlp solver

      // Auxiliaries should be correct. solution should be the one found
      // at the node even if not as good as the best known.

      // save violation flag and disregard it while adding cut at NLP
      // point (which are not violated by the current, NLP, solution)
      bool save_av = addviolated_;
      addviolated_ = false;

      problem_ -> domain () -> push 
	(problem_ -> nVars (), 
	 problem_ -> domain () -> x  (), 
	 problem_ -> domain () -> lb (), 
	 problem_ -> domain () -> ub ());
      // fill originals with nlp values
      CoinCopyN (nlpSol, problem_ -> nOrig (), problem_ -> domain () -> x ());
      problem_ -> initAuxs ();

      if (jnlst_ -> ProduceOutput (J_VECTOR, J_CONVEXIFYING)) {
	jnlst_ -> Printf(J_VECTOR, J_CONVEXIFYING,"== genrowcuts on NLP =============\n");
	for (int i = 0; i < problem_ -> nVars (); i++)
	  if (problem_ -> Var (i) -> Multiplicity () > 0)
	    jnlst_->Printf(J_VECTOR, J_CONVEXIFYING,"%4d %+20.8g [%+20.8g,%+20.8g]\n", i,
			   problem_ -> X  (i),
			   problem_ -> Lb (i),
			   problem_ -> Ub (i));
	jnlst_->Printf(J_VECTOR, J_CONVEXIFYING,"=============================\n");
      }

      genRowCuts (si, cs, nchanged, changed, //info, 
		  chg_bds, true);  // add cuts

      problem_ -> domain () -> pop (); // restore point

      addviolated_ = save_av;     // restore previous value

      //    if (!firstcall_) // keep solution if called from extractLinearRelaxation()
      babInfo -> setHasNlpSolution (false); // reset it after use //AW HERE
    } else {

      if (jnlst_ -> ProduceOutput (J_VECTOR, J_CONVEXIFYING)) {
	jnlst_ -> Printf(J_VECTOR, J_CONVEXIFYING,"== genrowcuts on LP =============\n");
	for (int i = 0; i < problem_ -> nVars (); i++)
	  if (problem_ -> Var (i) -> Multiplicity () > 0)
	    jnlst_->Printf(J_VECTOR, J_CONVEXIFYING,"%4d %+20.8g [%+20.8g,%+20.8g]\n", i,
			   problem_ -> X  (i),
			   problem_ -> Lb (i),
			   problem_ -> Ub (i));
	jnlst_->Printf(J_VECTOR, J_CONVEXIFYING,"=============================\n");
      }

      genRowCuts (si, cs, nchanged, changed, //info, 
		  chg_bds);
    }

    // change tightened bounds through OsiCuts
    if (nchanged)
      genColCuts (si, cs, nchanged, changed);

    if (firstcall_ && (cs.sizeRowCuts () >= 1))
      jnlst_->Printf(J_SUMMARY, J_CONVEXIFYING,
		     "Couenne: %d initial row cuts\n", cs.sizeRowCuts ());
  }

  // end of OBBT //////////////////////////////////////////////////////////////////////

  catch (int exception) {

    if ((exception == INFEASIBLE) && (!firstcall_)) {

      jnlst_ -> Printf (J_DETAILED, J_CONVEXIFYING,
			"generateCuts: Infeasible node detected\n");

      OsiColCut *infeascut = new OsiColCut;

      if (infeascut) {
	int i=0;
	double upper = -1., lower = +1.;
	infeascut -> setLbs (1, &i, &lower);
	infeascut -> setUbs (1, &i, &upper);
	cs.insert (infeascut);
	delete infeascut;
      }
    }

    if (babInfo) // set infeasibility to true in order to skip NLP heuristic
      babInfo -> setInfeasibleNode ();
  }

  delete [] chg_bds;
  free (changed);

  if (!firstcall_)
    problem_ -> domain () -> pop ();

  if (firstcall_) {

    fictitiousBound (cs, problem_, true);
    firstcall_  = false;
    ntotalcuts_ = nrootcuts_ = cs.sizeRowCuts ();
  }
  else ntotalcuts_ += (cs.sizeRowCuts () - nInitCuts);

  septime_ += CoinCpuTime () - now;

  if (jnlst_ -> ProduceOutput (J_VECTOR, J_CONVEXIFYING)) {

    if (cs.sizeColCuts ()) {
      jnlst_ -> Printf (J_DETAILED, J_CONVEXIFYING,"Couenne col cuts:\n");
      for (int i=0; i<cs.sizeColCuts (); i++) 
	cs.colCutPtr (i) -> print ();
    }

    jnlst_ -> Printf(J_VECTOR, J_CONVEXIFYING,"== on my way out of generateCuts =============\n");
    for (int i = 0; i < problem_ -> nVars (); i++)
      if (problem_ -> Var (i) -> Multiplicity () > 0)
	jnlst_->Printf(J_VECTOR, J_CONVEXIFYING,
		       "%4d %+20.8g [%+20.8g,%+20.8g] %+20.8g [%+20.8g,%+20.8g]\n", i,
		       problem_ -> X  (i),
		       problem_ -> Lb (i),
		       problem_ -> Ub (i),
		       si.getColSolution () [i],
		       si.getColLower    () [i],
		       si.getColUpper    () [i]);
    jnlst_->Printf(J_VECTOR, J_CONVEXIFYING,"=============================\n");
  }
}

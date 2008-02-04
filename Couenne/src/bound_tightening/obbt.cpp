/*
 * Name:    obbt.cpp
 * Author:  Pietro Belotti
 * Purpose: Optimality-Based Bound Tightening
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CglCutGenerator.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "CouenneSolverInterface.hpp"

#define OBBT_EPS 1e-3
#define MAX_OBBT_LP_ITERATION 100

// minimum #bound changed in obbt to generate further cuts
#define THRES_NBD_CHANGED 1

// maximum number of obbt iterations
#define MAX_OBBT_ITER 1

// defined in generateCuts.cpp
void sparse2dense (int ncols, t_chg_bounds *chg_bds, int *&changed, int &nchanged);


int CouenneProblem::call_iter (CouenneSolverInterface *csi, 
			       t_chg_bounds *chg_bds, 
			       const CoinWarmStart *warmstart, 
			       Bonmin::BabInfo *babInfo,
			       double *objcoe,
			       enum nodeType type,
			       int sense) const {

  int ncols   = csi -> getNumCols (),
      nimprov = 0;

  for (int ii=0; ii<ncols; ii++) {
    int i = evalOrder (ii);

    if ((Var (i) -> Type () == type) &&
	(Var (i) -> Multiplicity () > 0)) {

      int ni = obbt_iter (csi, chg_bds, warmstart, babInfo, objcoe, sense, i);

      if (ni < 0) return ni;
      nimprov += ni;
    }
  }

  return nimprov;
}


/// Optimality based bound tightening -- inner loop

int CouenneProblem::obbtInner (CouenneSolverInterface *csi,
			       OsiCuts &cs,
			       t_chg_bounds *chg_bds,
			       Bonmin::BabInfo * babInfo) const {

  // set large bounds to infinity (as per suggestion by JJF)

  int ncols = csi -> getNumCols ();
  const double *lb = csi -> getColLower (),
               *ub = csi -> getColUpper ();

  double inf = csi -> getInfinity ();

  for (int i=ncols; i--;) {
    if (lb [i] < - COUENNE_INFINITY) csi -> setColLower (i, -inf);
    if (ub [i] >   COUENNE_INFINITY) csi -> setColUpper (i,  inf);
  }

  //  csi -> setHintParam (OsiDoDualInResolve, false);

  // setup cloned interface for later use
  csi -> setObjSense (1); // minimization
  csi -> setIntParam (OsiMaxNumIteration, MAX_OBBT_LP_ITERATION);
  csi -> applyCuts (cs);   // apply all (row+column) cuts to formulation
  csi -> initialSolve ();

  const CoinWarmStart *warmstart = csi -> getWarmStart ();

  // improve each bound

  double *objcoe = (double *) malloc (ncols * sizeof (double));

  // set obj function coefficients to zero
  for (int i=ncols; i--;)
    *objcoe++ = 0.;
  objcoe -= ncols;

  csi -> setObjective (objcoe);
  csi -> setObjSense (1);        // minimization

  int nimprov = 0;
 
  const int Infeasible = 1;

  try {

    int ni;

    if ((ni = call_iter (csi, chg_bds, warmstart, babInfo, objcoe, VAR,  1)) < 0) throw Infeasible;
    nimprov += ni;

    if ((ni = call_iter (csi, chg_bds, warmstart, babInfo, objcoe, VAR, -1)) < 0) throw Infeasible;
    nimprov += ni;

    if ((ni = call_iter (csi, chg_bds, warmstart, babInfo, objcoe, AUX,  1)) < 0) throw Infeasible;
    nimprov += ni;

    if ((ni = call_iter (csi, chg_bds, warmstart, babInfo, objcoe, AUX, -1)) < 0) throw Infeasible;
    nimprov += ni;

    free (objcoe);
    delete warmstart;
  }

  catch (int exception) {
    if (exception == Infeasible)
      return -1;
  }

  return (nimprov);
}


// Optimality based bound tightening -- main loop
int CouenneProblem::obbt (const CouenneCutGenerator *cg,
			  const OsiSolverInterface &si,
			  OsiCuts &cs,
			  const CglTreeInfo &info,
			  Bonmin::BabInfo * babInfo,
			  int nchanged,
			  int *changed,
			  t_chg_bounds *chg_bds) {

  // Do OBBT if:

  if (doOBBT_ &&                        // flag is checked, AND
      ((logObbtLev_ != 0) ||               // (parameter is nonzero OR
       (info.level == 0)) &&               //  we are at root node), AND
      (info.pass == 0) &&               // at first round of cuts, AND 
      ((logObbtLev_ < 0) ||               // (logObbtLev = -1, OR
       (info.level <= logObbtLev_) ||     //  depth is lower than COU_OBBT_CUTOFF_LEVEL, OR
                                          //  probability inversely proportional to the level)
       (CoinDrand48 () < pow (2., (double) logObbtLev_ - (info.level + 1))))) {

    // TODO: why check info.pass==0? Why not more than one pass? It
    // should be anyway checked that info.level be >= 0 as <0 means
    // first call at root node

    CouenneSolverInterface *csi = dynamic_cast <CouenneSolverInterface *> (si.clone (true));

    csi -> setupForRepeatedUse ();

    int nImprov, nIter = 0;

    while ((nIter++ < MAX_OBBT_ITER) &&
	   ((nImprov = obbtInner (csi, cs, chg_bds, babInfo)) > 0)) {

      /// OBBT has tightened, add improved bounds
      sparse2dense (nVars(), chg_bds, changed, nchanged);
      cg -> genColCuts (*csi, cs, nchanged, changed);

      if (nImprov >= THRES_NBD_CHANGED) {

	// only generate new row cuts if improvents are enough
	int nCurCuts = cs.sizeRowCuts ();
	cg -> genRowCuts (*csi, cs, nchanged, changed, info, chg_bds);

	if (nCurCuts == cs.sizeRowCuts ())
	  break; // repeat only if new cuts available

      } else break;
    }

    delete csi;

    if (nImprov < 0)
      jnlst_->Printf(J_DETAILED, J_CONVEXIFYING,
		     "### infeasible node after OBBT\n");

    if (nImprov < 0)
      return -1;
  }

  return 0;
}

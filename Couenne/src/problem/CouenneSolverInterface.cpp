/*
 * Name:    CouenneSolverInterface.cpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: Implementation of the OsiSolverInterface::resolve () method 
 *
 * (C) Carnegie-Mellon University, 2006, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "OsiClpSolverInterface.hpp"
#include "CouenneProblem.hpp"
#include "CouenneSolverInterface.hpp"
#include "CglTreeInfo.hpp"

//#define DEBUG

CouenneSolverInterface::CouenneSolverInterface (CouenneCutGenerator *cg /*= NULL*/)
  :
  OsiClpSolverInterface(),
  cutgen_ (cg),
  knowInfeasible_(false),
  knowOptimal_(false)
{}

CouenneSolverInterface::CouenneSolverInterface (const CouenneSolverInterface &src)
  :
  OsiSolverInterface (src),
  OsiClpSolverInterface (src),
  cutgen_ (src.cutgen_),
  knowInfeasible_ (src.knowInfeasible_),
  knowOptimal_ (src.knowOptimal_) {}


/// Destructor
CouenneSolverInterface::~CouenneSolverInterface () {
  //  if (cutgen_)
  //    delete cutgen_;
}


/// Solve initial LP relaxation 
void CouenneSolverInterface::initialSolve () 
{
  knowInfeasible_ = false;
  knowOptimal_ = false;
  OsiClpSolverInterface::initialSolve ();
}

bool CouenneSolverInterface::isProvenPrimalInfeasible() const
{
  if (knowInfeasible_) {
    return true;
  }
  return OsiClpSolverInterface::isProvenPrimalInfeasible();
}

bool CouenneSolverInterface::isProvenOptimal() const
{
  if (knowOptimal_) {
    return true;
  }
  return OsiClpSolverInterface::isProvenOptimal();
}

/// Defined in Couenne/src/convex/generateCuts.cpp
void sparse2dense (int, t_chg_bounds *, int *&, int &);


/// Resolve an LP relaxation after problem modification
void CouenneSolverInterface::resolve () {
  // TODO: if NLP point available, add new cuts BEFORE resolving --
  // and decrease number of cutting plane iterations by one, to
  // balance it
  //printf("NumRows in resolve = %d\n", getNumRows());
  knowInfeasible_ = false;
  knowOptimal_ = false;
  OsiClpSolverInterface::resolve ();
  //printf("obj value in resolve = %e\n",getObjValue());
}

/// Create a hot start snapshot of the optimization process.
void CouenneSolverInterface::markHotStart () {
  //printf(">>>> markHotStart\n");
  // Using OsiClpSolverInterface doesn't work yet...
  OsiSolverInterface::markHotStart ();
}

/// Optimize starting from the hot start snapshot.
void CouenneSolverInterface::solveFromHotStart() {

  OsiSolverInterface::solveFromHotStart();

  knowInfeasible_ = false;
  knowOptimal_ = false;
  const int ncols = cutgen_ -> Problem () -> nVars ();

  cutgen_ -> Problem () -> update (getColSolution (),
				   getColLower    (),
				   getColUpper    ());

  // This vector contains variables whose bounds have changed due to
  // branching, reduced cost fixing, or bound tightening below. To be
  // used with malloc/realloc/free

  t_chg_bounds *chg_bds = new t_chg_bounds [ncols];

  OsiCuts cs;

  Bonmin::BabInfo *babInfo = dynamic_cast <Bonmin::BabInfo *> (getAuxiliaryInfo ());

  if (! (cutgen_ -> Problem () -> boundTightening (chg_bds, babInfo))) {

#ifdef DEBUG
    printf ("#### BT says infeasible before re-solve\n");
#endif
    knowInfeasible_ = true;
    return;
  }

  int *changed = NULL, nchanged;
  sparse2dense (ncols, chg_bds, changed, nchanged);

  // change tightened bounds through OsiCuts
  if (nchanged)
    cutgen_ -> genColCuts (*this, cs, nchanged, changed);

  const int nRowsBeforeRowCuts = getNumRows();
  //printf("NumRows before getRowCuts = %d\n", getNumRows());
  cutgen_ -> genRowCuts (*this, cs, nchanged, changed, CglTreeInfo(),
			 chg_bds, false);

  // Now go through the list of cuts and apply the column cuts
  // directly as changes on bounds
  while(cs.sizeColCuts()) {
    const OsiColCut& ccut = cs.colCut(0);
    const CoinPackedVector& lbs = ccut.lbs();
    int nele = lbs.getNumElements();
    const int* idxs = lbs.getIndices();
    const double* eles = lbs.getElements();
    const double* bnds = getColLower();
    for (int i=0; i<nele; i++) {
      if (bnds[*idxs] < *eles) {
	setColLower(*idxs,*eles);
      }
      idxs++;
      eles++;
    }
    const CoinPackedVector& ubs = ccut.ubs();
    nele = ubs.getNumElements();
    idxs = ubs.getIndices();
    eles = ubs.getElements();
    bnds = getColUpper();
    for (int i=0; i<nele; i++) {
      if (bnds[*idxs] > *eles) {
	setColUpper(*idxs,*eles);
      }
      idxs++;
      eles++;
    }
    cs.eraseColCut(0); 
  }

  applyCuts (cs);

  const int nRowsAfterRowCuts = getNumRows();
  //printf("NumRows after applyCuts = %d\n", getNumRows());

  resolve();
  if (isProvenPrimalInfeasible()) {
    knowInfeasible_ = true;
  }
  if (isProvenOptimal()) {
    knowOptimal_ = true;
  }
  //printf("obj value = %e\n",getObjValue());

  // now undo the row cuts
  int nrowsdel = nRowsAfterRowCuts-nRowsBeforeRowCuts;
  int* rowsdel = new int[nrowsdel];
  for(int i=0; i<nrowsdel; i++) {
    rowsdel[i] = nRowsBeforeRowCuts+i;
  }
  deleteRows(nrowsdel, rowsdel);
  delete [] rowsdel;
  //printf("NumRows after deleting = %d\n", getNumRows());
}

/// Delete the hot start snapshot.
void CouenneSolverInterface::unmarkHotStart() {
  //printf("<<<< unmarkHotStart\n");
  OsiSolverInterface::unmarkHotStart();
}

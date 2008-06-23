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
  knowOptimal_ (src.knowOptimal_) 
{}

/// Destructor
CouenneSolverInterface::~CouenneSolverInterface () {
  //  if (cutgen_)
  //    delete cutgen_;
}


/// Solve initial LP relaxation 
void CouenneSolverInterface::initialSolve () {
  /*printf ("------------------------------------- INITIAL SOLVE\n");
  for (int i=0; i<getNumCols(); i++)
    printf ("%4d. %20.5g [%20.5g %20.5g]\n", 
	    i, getColSolution () [i],
	    getColLower () [i],
	    getColUpper () [i]);

	    cutgen_ -> Problem () -> print ();*/

  knowInfeasible_ = false;
  knowOptimal_ = false;

  OsiClpSolverInterface::initialSolve ();
  //writeLp ("initialLP");
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
  /*printf ("------------------------------------- RESOLVE\n");
  for (int i=0; i<getNumCols(); i++) 
    printf ("%4d. %20.5g [%20.5g %20.5g]\n", 
	    i, getColSolution () [i],
	    getColLower () [i],
	    getColUpper () [i]);*/

  // CUT THIS! Some about-to-be-resolved problems have variables with
  // lower +inf. That is, between CbcModel::initialSolve() and
  // CbcModel::resolve(). I couldn't spot where in Couenne this
  // happens.

  ////////////////////////////////////// Cut {
  /*const double 
    *lb = getColLower (),
    *ub = getColUpper ();

  for (int i=getNumCols(); i--;) {
    if (lb [i] >  COUENNE_INFINITY)
      setColLower (i, cutgen_ -> Problem () -> Lb (i));
    //setColLower (i, -COIN_DBL_MAX);//cutgen_ -> Problem () -> Lb (i));
    if (ub [i] < -COUENNE_INFINITY)
      setColUpper (i, cutgen_ -> Problem () -> Ub (i));
    //setColUpper (i,  COIN_DBL_MAX);//cutgen_ -> Problem () -> Ub (i));
    }*/
  ////////////////////////////////////// Cut }

  // TODO: if NLP point available, add new cuts BEFORE resolving --
  // and decrease number of cutting plane iterations by one, to
  // balance it
  //printf("NumRows in resolve = %d\n", getNumRows());

  //static int count = 0;
  //char filename [30];
  //sprintf (filename, "presol_%d", count);
  //writeLp (filename);

  //printf ("----------------------------- count = %d [%s]\n", count, filename);

  knowInfeasible_ = false;
  knowOptimal_    = false;

  OsiClpSolverInterface::resolve ();

  //sprintf (filename, "postsol_%d", count++);
  //writeLp (filename);

  /*printf("obj value in resolve = %e\n",getObjValue());
  printf ("after resolve, %p --> [%g,%g]\n", this,
	  getColLower () [getNumCols () - 1],
	  getColUpper () [getNumCols () - 1]);*/
}


/// Create a hot start snapshot of the optimization process.
void CouenneSolverInterface::markHotStart () {
  //printf(">>>> markHotStart\n");
  // Using OsiClpSolverInterface doesn't work yet...
  OsiSolverInterface::markHotStart ();
}


/// Delete the hot start snapshot.
void CouenneSolverInterface::unmarkHotStart() {
  //printf("<<<< unmarkHotStart\n");
  OsiSolverInterface::unmarkHotStart();
}



/// Optimize starting from the hot start snapshot.
void CouenneSolverInterface::solveFromHotStart() {

  //OsiClpSolverInterface::solveFromHotStart ();

  //#if 0
  knowInfeasible_ = false;
  knowOptimal_ = false;

  /*
  const int ncols = cutgen_ -> Problem () -> nVars ();

  cutgen_ -> Problem () -> domain () -> push
    (cutgen_ -> Problem () -> nVars (),
     getColSolution (),
     getColLower    (),
     getColUpper    ());

  // This vector contains variables whose bounds have changed due to
  // branching, reduced cost fixing, or bound tightening below. To be
  // used with malloc/realloc/free

  t_chg_bounds *chg_bds = new t_chg_bounds [ncols];

  OsiCuts cs;

  Bonmin::BabInfo *babInfo = dynamic_cast <Bonmin::BabInfo *> (getAuxiliaryInfo ());

  if (cutgen_ -> Problem () -> doFBBT () && 
      (!(cutgen_ -> Problem () -> boundTightening (chg_bds, babInfo)))) {

#ifdef DEBUG
    printf ("#### BT says infeasible before re-solve\n");
#endif
    knowInfeasible_ = true;
    cutgen_ -> Problem () -> domain () -> pop ();
    return;
  }

  int *changed = NULL, nchanged;
  sparse2dense (ncols, chg_bds, changed, nchanged);

  // change tightened bounds through OsiCuts
  if (nchanged)
    cutgen_ -> genColCuts (*this, cs, nchanged, changed);

  const int nRowsBeforeRowCuts = getNumRows();
  //printf("NumRows before getRowCuts = %d\n", getNumRows());
  cutgen_ -> genRowCuts (*this, cs, nchanged, changed, //CglTreeInfo(),
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
	//printf ("setcolLower %d %g\n", *idxs, *eles);
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
	//printf ("setcolUpper %d %g\n", *idxs, *eles);
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
  */

  resolve();

  if (isProvenPrimalInfeasible ()) knowInfeasible_ = true;
  if (isProvenOptimal ())          knowOptimal_    = true;

  //printf("obj value = %e\n",getObjValue());

  // now undo the row cuts
  /*
  int nrowsdel = nRowsAfterRowCuts-nRowsBeforeRowCuts;
  int* rowsdel = new int[nrowsdel];
  for(int i=0; i<nrowsdel; i++) {
    rowsdel[i] = nRowsBeforeRowCuts+i;
  }
  deleteRows(nrowsdel, rowsdel);
  delete [] rowsdel;
  //printf("NumRows after deleting = %d\n", getNumRows());

  cutgen_ -> Problem () -> domain () -> pop ();
  */
  //#endif
}

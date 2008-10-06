/*
 * Name:    OsiCuts2MatrVec.cpp
 * Author:  Pietro Belotti
 * Purpose: turn OsiCuts objects into coefficient matrix and rhs vector
 *
 * (C) Carnegie-Mellon University, 2008. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneDisjCuts.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinHelperFunctions.hpp"
#include "CouennePrecisions.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"


// add CGLP columns to solver interface; return number of columns
// added (for later removal)
int CouenneDisjCuts::OsiCuts2MatrVec (OsiSolverInterface *cglp,
				      OsiCuts *cuts,
				      int displRow,
				      int displRhs) const {
  int
    ncols  = 0,
    nrcuts = cuts -> sizeRowCuts (),
    nccuts = cuts -> sizeColCuts (),
    ncgrow = cglp -> getNumRows () - 1,
    nnzR   = 0,
    ncC    = 0;

  if (!(nrcuts || nccuts))
    return 0;

  // count nonzero in row cuts
  for (int i=nrcuts; i--;) {

    OsiRowCut *cut = cuts -> rowCutPtr (i);

    if ((cut -> sense () == 'E') ||
	(cut -> sense () == 'R')) {

      nnzR += 2 * cut -> row (). getNumElements ();
      nrcuts++; // can increase as not used

    } else nnzR += cut -> row (). getNumElements ();
  }

  // count bound constraints in column cuts
  for (int i=nccuts; i--;) {
    OsiColCut *cut = cuts -> colCutPtr (i);
    ncC += 
      cut -> lbs ().getNumElements () + 
      cut -> ubs ().getNumElements ();
  }

  int 
     nnz     = 2 * (nnzR + 2*nrcuts + 3*ncC),
    *indices = new int [nnz],
    *start   = new int [nrcuts + ncC + 1],
     curel   = 0;

  double 
    *elements = new double [nnz],               // for row cuts + col cuts
    *collb    = new double [2*(nrcuts + ncC)],  // lower bounds for new columns
    *colub    = new double [2*(nrcuts + ncC)],  // upper 
    *obj      = new double [2*(nrcuts + ncC)];  // objective coefficient (zero)

  // trivial, lower/upper bounds and objective coefficients
  CoinFillN (collb, 2*(nrcuts + ncC), 0.);
  CoinFillN (colub, 2*(nrcuts + ncC), 1.);
  CoinFillN (obj,   2*(nrcuts + ncC), 0.);

  // scan OsiColCuts ////////////////////////////////////////

  double *saveEl  = elements;
  int    *saveInd = indices;

  for (int i = nccuts; i--;) {

    OsiColCut *cut = cuts -> colCutPtr (i);

    const CoinPackedVector
      &lbs = cut -> lbs (),
      &ubs = cut -> ubs ();

    // lower bounds
    const int
      *lind = lbs. getIndices (), 
       nlcc = lbs. getNumElements ();
    const double *lele = lbs. getElements ();

    for (int j = nlcc; j--; lind++, lele++)

      if (couenneCG_ -> Problem () -> Var (*lind) -> Multiplicity () > 0) {
	*start++ = curel;
	*elements++ = -1.;        *indices++ = displRow + *lind;
	if (fabs (*lele) > COUENNE_EPS) 
	  {*elements++ = -*lele;     *indices++ = displRhs; curel++;}
	*elements++ =  1.;        *indices++ = ncgrow;
	curel += 2;
      }

    // upper bounds
    const int
      *uind = ubs. getIndices (), 
       nucc = ubs. getNumElements ();
    const double *uele = ubs. getElements ();

    for (int j = nucc; j--; uind++, uele++)
      if (couenneCG_ -> Problem () -> Var (*uind) -> Multiplicity () > 0) {
	*start++ = curel;
	*elements++ =  1.;        *indices++ = displRow + *uind;
	if (fabs (*uele) > COUENNE_EPS) 
	  {*elements++ = *uele;     *indices++ = displRhs; curel++;}
	*elements++ =  1.;        *indices++ = ncgrow;
	curel += 2;
      }

    ncols += nlcc + nucc;
  }

  elements = saveEl;
  indices  = saveInd;
    
  //elements -= (3 * ncols);
  //indices  -= (3 * ncols);
  start    -= ncols;

  start [ncols] = curel; // may go


  if (jnlst_ -> ProduceOutput (J_MATRIX, J_DISJCUTS)) {
    printf ("%d cuts, have %d cols and cur el is %d. Now for the %d row cuts\n",
	    nccuts, ncols, curel, nrcuts);

    printf ("matrix (.. %d) %d elements:\n", ncols, curel);

    printf ("start: "); for (int i=0; i<=ncols; i++) printf ("%d ", start [i]);

    printf ("\nElements:\n"); 
    for (int i=0, j=0; i<ncols; i++) {
      for (int k=0; k<start[i+1] - start[i]; k++, j++) 
	printf ("(%d %g) ", indices [j], elements [j]);
      printf ("\n");
    }
  }

  // scan OsiRowCuts /////////////////////////////////////////////

  for (int i = cuts -> sizeRowCuts (); i--;) {

    OsiRowCut *cut = cuts -> rowCutPtr (i);

    const CoinPackedVector &row = cut -> row ();

    const double 
      rhs    = cut -> rhs   (),
      rng    = cut -> range (),
      *rowEl = row. getElements ();

    const int 
      *rowIn = row. getIndices (),
      rowNE  = row. getNumElements ();

    switch (cut -> sense ()) {

    case 'L': 

      start [ncols++] = curel;
      CoinCopyDisp (rowIn, rowNE, indices  + curel, displRow);
      CoinCopyN    (rowEl, rowNE, elements + curel);
      curel += rowNE;
      if (fabs (rhs) > COUENNE_EPS) {indices  [curel] = displRhs;  elements [curel++] = rhs;}
      indices  [curel] = ncgrow;    elements [curel++] = 1.;

      break;

    case 'E':
    case 'R':

      start [ncols++] = curel;
      CoinCopyDisp (rowIn, rowNE, indices  + curel, displRow);
      CoinCopyN    (rowEl, rowNE, elements + curel);
      curel += rowNE;
      if (fabs (rhs+rng) > COUENNE_EPS) {indices  [curel] = displRhs; elements [curel++] = rhs + rng;}
      // rng only used here, zero for 'E'
      indices  [curel] = ncgrow;   elements [curel++] = 1.;

      start [ncols++] = curel;
      CoinCopyDisp (rowIn, rowNE, indices  + curel, displRow);
      CoinInvN     (rowEl, rowNE, elements + curel);
      curel += rowNE;
      if (fabs (rhs) > COUENNE_EPS) {indices  [curel] = displRhs;  elements [curel++] = -rhs;}
      indices  [curel] = ncgrow;    elements [curel++] = 1.;

      break;

    case 'G':
      start [ncols++] = curel;
      CoinCopyDisp (rowIn, rowNE, indices  + curel, displRow);
      CoinInvN     (rowEl, rowNE, elements + curel);
      curel += rowNE;
      if (fabs (rhs) > COUENNE_EPS) {indices  [curel] = displRhs;  elements [curel++] = -rhs;}
      indices  [curel] = ncgrow;    elements [curel++] = 1.;

      break;

    default: printf ("unknown type of cut\n");
      exit (-1);
    }
  }

  start [ncols] = curel;

  //  ncols=1;start[1]=3;
  /*  *start=0; start[1]=1;
  *indices=0;
  *elements=4.55;
  *collb=-3;
  *colub=+3;
  *obj=4;*/


  if (jnlst_ -> ProduceOutput (J_MATRIX, J_DISJCUTS)) {
    printf ("===================\nmatrix (.. %d) %d elements:\n", ncols, curel);

    printf ("start: "); for (int i=0; i<=ncols; i++) printf ("%d ", start [i]);

    printf ("\nElements:\n"); 
    for (int i=0, j=0; i<ncols; i++) {
      for (int k=0; k<start[i+1] - start[i]; k++, j++) 
	printf ("(%d %g) ", indices [j], elements [j]);
      printf ("\n");
    }
  }

  /*{
    const CoinPackedMatrix *m = cglp->getMatrixByCol();
    printf ("before: size_ = %d, start [%d] = %d\n", 
	    m -> getNumElements (), 
	    m -> getSizeVectorLengths(),
	    m -> getVectorStarts () [m -> getSizeVectorLengths()]);
	    }*/

  cglp -> addCols (ncols,    start,
		   indices,  elements,
		   collb,    colub,   
		   obj);

  /*{
    const CoinPackedMatrix *m = cglp->getMatrixByCol();
    printf ("after: size_ = %d, start [%d] = %d\n", 
	    m -> getNumElements (), 
	    m -> getSizeVectorLengths(),
	    m -> getVectorStarts () [m -> getSizeVectorLengths()]);
	    }*/

  delete [] elements;
  delete [] collb;
  delete [] colub;
  delete [] obj;

  delete [] indices;
  delete [] start;

  return ncols;
}

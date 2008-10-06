/*
 * Name:    disjCut.cpp
 * Author:  Pietro Belotti
 * Purpose: generate one disjunctive cut based on a single disjunction
 *
 * (C) Carnegie-Mellon University, 2008. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinTime.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CouennePrecisions.hpp"
#include "CouenneDisjCuts.hpp"

#define COEFF_BOUNDS 1.e10

#define MIN_NUM_COEFF 1.e-9
#define MAX_NUM_COEFF 1.e+9
#define MAX_NUM_RATIO 1.e+6


// print content of sparse matrix
void printMatrix (int nrows, int ncols, int nel, 
		  const int *start, const int *len, 
		  const int *ind, const double *el);

// same with CoinPackedMatrix
void printMatrix   (const CoinPackedMatrix *A);

// same with si.GetMatrixByRow()
void printLPMatrix (const OsiSolverInterface &si);

// add columns specified by start/len/ind/el to matrix Astd
void addSubMatr (int *start, int *len, int *ind, double *el, 
		 CoinPackedMatrix &Astd, CoinPackedVector &rstd, 
		 int &cur, int &curCol, int dispM, int dispVec, int nrows);


/// generate one disjunctive cut from one CGLP
int CouenneDisjCuts::generateDisjCuts (std::vector <std::pair <OsiCuts *, OsiCuts *> > &disjunctions, 
				       OsiSolverInterface &si, 
				       OsiCuts &cs, 
				       const CglTreeInfo &info) const {

  // create CGLP with si+left and si+right, for each (left,right) in
  // the vector of disjunctions
  //
  // Here's the CGLP (row and column order as shown)
  //
  //        {n}          {1}    {m}    {m}   {mL} {mR}
  //
  // max    alpha xbar - beta                                        --- maximize so can copy xbar
  // s.t.  -alpha              + u A       + u'C           =  0      --- n rows
  //       -alpha                    + v A       + v'D     =  0      --- n
  //                    -beta  + u b       + u'c          <=  0      --- 1
  //                    -beta        + v b       + v'd    <=  0      --- 1
  //                           |(u,    v,    u',   v')|_1  =  1      --- normalized multipliers
  //
  //                             u,    v,    u',   v'     >=  0      --- non-negativity 
  //
  //
  //  And here are the submatrices (M' is M transposed)
  //
  //     First        Second  Third  Fourth Fifth Sixth   --- Columns
  //
  //       xbar          -1       0     0     0     0
  //      ---------------------------------------------
  //       -I_n           .       A'    .     C'    .     = 0
  //       -I_n           .       .     A'    .     D'    = 0
  //        .            -1       b'    .     c'    .    <= 0
  //        .            -1       .     b'    .     d'   <= 0
  //        .             .       e     e     e     e     = 1
  //
  //
  // build a different A such that
  // 
  // - only active rows (resp. columns) of A are included if active_rows_
  //   (resp. active_columns_) are set
  // - equality constraints are duplicated in a <= and a >= constraint
  // - >= constraints are inverted
  //
  // Also, add disjunctive cut to CGLP for use with next disjunction

  // put matrix from base problem in canonical form //////////////////////////////
  CoinPackedMatrix Astd;
  CoinPackedVector rstd;
  OsiSI2MatrVec (Astd,  rstd,  si); 

  int
    n   = si.   getNumCols (),      //mC   = 2*n + 3,
    m   = Astd. getMajorDim (),     nC   = 1 + n + 2 * m,
    nnz = Astd. getNumElements (),  nnzC = 2 * (n + 1 + nnz + 2*m);

  if (jnlst_ -> ProduceOutput (J_DETAILED, J_DISJCUTS))
    printf ("canonical form has %d cols, %d rows, %d nonzeros --> cglp has %d,%d,%d\n", 
	    n, m, nnz, nC, 2*n + 3, nnzC);

  double 
    *elements = new double [nnzC];

  int 
    *indices = new int [nnzC],
    *start   = new int [nC + 1],
    *length  = new int [nC],
    cur      = 0,
    curCol   = 0;

  // first column: two identity matrices
  for (int i=0, i2 = n; i<n;) {
    start   [curCol]   = cur;
    length  [curCol++] = 2;
    indices [cur] =    i++; elements [cur++] = -1.;
    indices [cur] =   i2++; elements [cur++] = -1.;
  }

  // second column: two "-1" at position 2n and 2n+1
  start   [curCol]   = cur;
  length  [curCol++] = 2;
  indices [cur] = 2*n;   elements [cur++] = -1.;
  indices [cur] = 2*n+1; elements [cur++] = -1.;

  // third...
  addSubMatr (start + curCol, length + curCol, 
	      indices + cur, elements + cur, 
	      Astd, rstd,
	      cur, curCol, 
	      0, 2*n,   2*n+2);

  if (jnlst_ -> ProduceOutput (J_MATRIX, J_DISJCUTS)) {
    printf ("with third column\n");
    printMatrix (curCol, 2*n+3, cur, start, length, indices, elements);
  }

  // ... and fourth column: get single rows from Astd
  addSubMatr (start + curCol, length + curCol, 
	      indices + cur, elements + cur, 
	      Astd, rstd, 
	      cur, curCol, 
	      n, 2*n+1, 2*n+2);

  if (jnlst_ -> ProduceOutput (J_MATRIX, J_DISJCUTS)) {
    printf ("with 4th column\n");
    printMatrix (curCol, 2*n+3, cur, start, length, indices, elements);
  }

  CoinPackedMatrix *baseA = new CoinPackedMatrix;

  baseA -> assignMatrix (true,     // column ordered
			 2*n+3,    // minor dimension
			 curCol,   // major dimension
			 cur,      // number of elements
			 elements, // elements
			 indices,  // indices
			 start,    // starting positions
			 length);  // length

  //printf ("should be copy of above\n");
  //printMatrix (baseA);

  OsiClpSolverInterface cglp;

  cglp. messageHandler () -> setLogLevel (0);

  int 
    N = baseA -> getMajorDim (),        // # cols in base problem
    M = baseA -> getMinorDim ();        // # rows in base problem

  assert (M == 2 * n + 3);

  // vectors of the problem
  double
    *collb  = new double [N], // variable lower bounds
    *colub  = new double [N], // variable upper bounds
    *obj    = new double [N], // objective coefficients
    *rowrhs = new double [M], // right hand sides (all zero except the last, 1)
    *rowrng = new double [M]; // row range (empty)

  // bounds
  CoinFillN (collb,       n+1,       -COEFF_BOUNDS);
  CoinFillN (collb + n+1, N - (n+1),  0.);
  CoinFillN (colub,       n+1,        COEFF_BOUNDS);
  CoinFillN (colub + n+1, N - (n+1),  1.);

  // objective coefficients
  CoinCopyN (si.getColSolution (), n, obj);
  obj [n] = -1.;
  CoinFillN (obj + (n+1), N-(n+1),    0.);

  // rhs
  CoinFillN (rowrhs,      M-1,        0.);
  rowrhs [M-1] = 1.;

  // rhs range
  CoinFillN (rowrng,      M,          COIN_DBL_MAX);

  // signs of the inequalities
  char *rowsen = new char [M];
  CoinFillN (rowsen, M, 'E');
  rowsen [M-3] = rowsen [M-2] = 'L';

  cglp.assignProblem (baseA,   // matrix
		      collb,   // lower bounds 
		      colub,   // upper bounds
		      obj,     // obj coefficients
		      rowsen,  // row sense
		      rowrhs,  // right hand sides
		      rowrng); // no row range

  // this is a maximization problem
  cglp. setObjSense (-1);

  /////////////////////////////////////////////////////////////////

  // generate and solve one CGLP for each disjunction

  bool first = true;

  for (std::vector <std::pair <OsiCuts *, OsiCuts *> >::iterator disjI = disjunctions.begin ();
       (disjI != disjunctions.end ()) && (CoinCpuTime () < cpuTime_); ++disjI) {

    OsiCuts
      *left  = disjI -> first,
      *right = disjI -> second;

    int 
      ncolLeft  = OsiCuts2MatrVec (&cglp,  left, 0, 2*n),
      ncolRight = OsiCuts2MatrVec (&cglp, right, n, 2*n+1);    

    /*char filename [30];
    static int iter = 0;
    sprintf (filename, "cglp-%04d-%04d-%04d", info.level, iter++, info.pass);
    cglp.writeLp (filename);*/

    if (jnlst_ -> ProduceOutput (J_MOREMATRIX, J_DISJCUTS)) {
      printf ("current CGLP:\n");
      printLPMatrix (cglp);
    }

    if (first) {cglp.initialSolve (); first = false;}
    else        cglp.resolve ();

    if (cglp. isProvenOptimal () && (cglp.getObjValue () > COUENNE_EPS)) {

      const double *AlphaBeta = cglp. getColSolution ();

      int    *colInd = NULL, nnzCut = 0;
      double *colCoe = NULL;

      // count nonzero entries, compute ratio max/min coefficient
      double mincoeff = COIN_DBL_MAX, maxcoeff = 0.;

      for (register int i=n+1; i--;) {
	double value = fabs (AlphaBeta [i]);
	if (value == 0.) continue;
	if (value > maxcoeff) maxcoeff = value;
	if (value < mincoeff) mincoeff = value;
	if ((maxcoeff            > MAX_NUM_COEFF) ||
	    (maxcoeff            < MIN_NUM_COEFF) ||
	    (maxcoeff / mincoeff > MAX_NUM_RATIO)) 
	  break;
	nnzCut++;
      }

      if (nnzCut &&
	  (maxcoeff            < MAX_NUM_COEFF) &&
	  (maxcoeff            > MIN_NUM_COEFF) &&
	  (maxcoeff / mincoeff < MAX_NUM_RATIO)) {

	// cut data
	double *nzcoeff = new double [nnzCut];
	int    *indices = new int    [nnzCut];

	// fill in indices and coefficient
	for (int i = nnzCut = 0; i<n; i++)
	  if (fabs (AlphaBeta [i]) > MIN_NUM_COEFF) {
	    indices [nnzCut]   = i;
	    nzcoeff [nnzCut++] = AlphaBeta [i];
	  }

	OsiRowCut *cut = new OsiRowCut;
	cut -> setRow (nnzCut, indices, nzcoeff);
	cut -> setUb  (AlphaBeta [n]);

	/*if (1) {

	  printf ("---- RESOLVING\n");
	  si.applyRowCuts (1, cut);
	  si.writeLp ("added");
	  si.resolve ();
	  printf ("---- RESOLVED\n");

	  double *obj    = new double [N]; // objective coefficients

	  // objective coefficients
	  CoinCopyN (si.getColSolution (), n, obj);
	  obj [n] = -1.;
	  CoinFillN (obj + (n+1), N-(n+1),    0.);

	  cglp.setObjective (obj);
	  }*/
 
	// add it to CGLP
	if (addPreviousCut_) {

	  colInd = new int    [2 * (nnzCut + 2)];
	  colCoe = new double [2 * (nnzCut + 2)];

	  // first column
	  CoinCopyN    (nzcoeff, nnzCut, colCoe);
	  CoinCopyN    (indices, nnzCut, colInd); 
	  colInd [nnzCut]       = 2*n;   colCoe [nnzCut]   = AlphaBeta [n];
	  colInd [nnzCut+1]     = 2*n+2; colCoe [nnzCut+1] = 1; // entry in norm constraint

	  // second column
	  CoinCopyN    (nzcoeff, nnzCut, colCoe + nnzCut + 2);
	  CoinCopyDisp (indices, nnzCut, colInd + nnzCut + 2, n); 
	  colInd [2*nnzCut + 2] = 2*n+1; colCoe [2*nnzCut+2] = AlphaBeta [n];
	  colInd [2*nnzCut + 3] = 2*n+2; colCoe [2*nnzCut+3] = 1.; // entry in norm constraint

	  // extra vectors
	  double lb  [2] = {0., 0.};
	  double ub  [2] = {1., 1.};
	  double obj [2] = {0., 0.};

	  int start [3];
	  *start = 0;
	  start [2] = 2 * (start [1] = nnzCut + 2);

	  cglp. addCols (2,        // const int numcols, 
			 start,    // const int* columnStarts,
			 colInd,   // const int* rows, 
			 colCoe,   // const double* elements,
			 lb,       // const double* collb, 
			 ub,       // const double* colub,   
			 obj);     // const double* obj

	  delete [] colCoe;
	  delete [] colInd;
	}

	delete [] nzcoeff;
	delete [] indices;

	if (jnlst_ -> ProduceOutput (J_DETAILED, J_DISJCUTS)) {
	  printf ("====== disjunctive cut: "); 
	  cut -> print ();
	}

	// add cut to cs
	cs. insert (cut);
      }
    }

    // remove last ncolLeft + ncolRight columns from cglp
    int *delIndices = new int [ncolLeft + ncolRight];
    for (register int nc = ncolLeft + ncolRight, j = N + nc; nc--;)
      *delIndices++ = --j;
    delIndices -= (ncolLeft + ncolRight);
    cglp.deleteCols (ncolLeft + ncolRight, delIndices);
    delete [] delIndices;
  }

  return COUENNE_FEASIBLE;
}

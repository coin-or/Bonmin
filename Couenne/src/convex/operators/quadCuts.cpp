/*
 * Name:    quadCuts.cpp
 * Author:  Pierre Bonami
 * Purpose: based on upper and lower convexification, add cuts to convexify
 *
 * (C) International Business Machines 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "exprQuad.hpp"

#include "CouenneProblem.hpp"
#include "CouenneCutGenerator.hpp"

#include "CoinHelperFunctions.hpp"

//#define DEBUG

void exprQuad::quadCuts (exprAux *w, OsiCuts &cs, const CouenneCutGenerator *cg){

  assert (dIndex_ != NULL);

#ifdef DEBUG
  std::cout<<"Expression has "<<nlterms_<<" linear terms and "
           <<nqterms_<<" quadratic terms. " << std::endl;

  printf ("Q\n");
  for (int i=0; i<nqterms_; i++)
    printf ("<%d,%d,%g>\n",  qindexI_ [i], qindexJ_ [i], qcoeff_ [i]);

  printf ("b\n");
  for (int i=0; i < nlterms_; i++)
    printf ("<%d,%g>\n",  index_ [i], coeff_ [i]);

  if (c0_) 
    printf ("<c0 = %g>\n", c0_);

  if (dCoeffLo_ && dCoeffUp_ && dIndex_) {
    printf ("alpha\n");
    for (int i=0; i<nDiag_; i++)
      printf ("[%d,%g,%g]\n", dIndex_ [i], dCoeffLo_ [i], dCoeffUp_ [i]);
  }
#endif

  // TODO: compute (*this') () where this' is convexification. If
  // current point is between nonconvex quad form and its
  // convexification, a cut won't help

  // Get on which side constraint is violated to get the good lambda

  double 
     varVal  = (*w)    (), 
     exprVal = (*this) (),
    *lambda  = (varVal < exprVal) ? 
      dCoeffLo_ : // Use under-estimator
      dCoeffUp_,  // Use  over-estimator
     convVal = 0;

  const CouenneProblem& problem = *(cg -> Problem ());
  const int numcols = problem.nVars ();

  const double 
    *colsol = problem.X  (), // current solution
    *lower  = problem.Lb (), //         lower bound
    *upper  = problem.Ub (); //         upper

  // compute lower or upper convexification and check if it contains
  // the current point

  if (dIndex_) {

    convVal = exprVal;

    if (lambda) { // there is a convexification, check if out of
		  // current point

      for (int i=0; i<nDiag_; i++) {

	int ind = dIndex_ [i];
	CouNumber xi = colsol [ind];
	convVal += lambda [i] * (xi - lower [ind]) * (upper [ind] - xi);
      }

      if (varVal < exprVal) {if (convVal < varVal) return;}
      else                  {if (convVal > varVal) return;}
    }
  }

#ifdef DEBUG
  std::cout << "Point to cut: "; for (int i = 0 ; i < numcols ; i++) std::cout << colsol [i] << ", ";
  printf (" (w,f(x),c) = (%g,%g,%g) -- lambda = [", (*w) (), exprVal, convVal);
  for (int i=0; i < nDiag_ ; i++) printf ("%g ", lambda [i]); printf ("]\n");
#endif

  // Initialize by copying $a$ into a dense vector and computing Q x^*
  double * Qxs = new double [numcols]; // sparse coefficient vector, $Qx^*$
  CoinFillN (Qxs, numcols, 0.);

  // Compute 2 * Q x^*.
  for (int k = 0 ; k < nqterms_ ; k++) {

    int qi = qindexI_ [k],
        qj = qindexJ_ [k];

    CouNumber qc = qcoeff_ [k];

    if (qi != qj) {
      Qxs [qi] += 2 * qc * colsol [qj]; // contribution of element $q_{ij}$ to (Qx)_i
      Qxs [qj] += 2 * qc * colsol [qi]; //                         $q_{ij}$    (Qx)_j
    }
    // elements on the diagonal are not halved upon reading
    else Qxs [qi] += 2 * qc * colsol [qi];
  }

#ifdef DEBUG
  printf ("2Qx = ("); for (int i = 0; i < numcols; i++) printf ("%g ", Qxs [i]); printf (")\n");
#endif

  // Add a to it.
  for (int i = 0 ; i < nlterms_ ; i++)
    Qxs [index_ [i]] += coeff_ [i];

  // multiply Qx^* by x^*^T again and store the result for the lower
  // bound into constant term

  double a0 = 0;//- c0_; // constant term

  for (int i = 0; i < numcols; i++){
    a0 += Qxs [i] * colsol [i];
    //    Qxs [i] *= 2;
  }

  // And myself
  Qxs [w -> Index ()] -= 1;

#ifdef DEBUG
  printf ("2Qx = ("); for(int i = 0; i < numcols; i++) printf ("%g ", Qxs [i]); printf (")[%g]\n",a0);
#endif

  a0 -= exprVal;

  if (lambda != NULL) // Now the part which depends on lambda, if there is one

    for (int k = 0 ; k < nDiag_ ; k++) {

      int ind = dIndex_ [k];
      CouNumber xi = colsol [ind],
	coeff = lambda [k] * (lower  [ind] + upper  [ind] - 2 * xi);

      a0 += coeff * xi - lambda [k] * (xi - lower [ind]) * (upper [ind] - xi);
      //      a0 += lambda [k] * lower  [ind] * upper  [ind];
      //      a0 -= lambda [k] * colsol [ind] * colsol [ind];

      Qxs [ind] += coeff;
      //Qxs [ind] -= lambda [k] * (colsol [ind]) * 2;
    }

  // Count the number of non-zeroes
  int nnz = 0;
  for (int i=0; i < numcols ; i++)
    if (fabs (Qxs [i]) > COUENNE_EPS)
      nnz++;

#ifdef DEBUG
  printf ("2Qx = (");
  for(int i = 0; i < numcols; i++)
    printf ("%g ", Qxs [i]);
  printf (")[%g], %d nz\n", a0, nnz);
#endif

  // Pack the vector into a CoinPackedVector and generate the cut.
  CoinPackedVector a (false);
  a.reserve (nnz);

#ifdef DEBUG
  CouNumber lhs = 0, lhsc = 0,
    *optimum = cg -> Problem () -> bestSol (),
    *current = cg -> Problem () -> X ();
#endif

  for (int i = 0 ; i < numcols ; i++)

    if (fabs (Qxs [i]) > COUENNE_EPS) {

      // compute violation
#ifdef DEBUG
      if (optimum) {
	printf ("%+g * %g ", Qxs [i], optimum [i]);
	lhs  += Qxs [i] * optimum [i];
      }
      lhsc += Qxs [i] * current [i];
#endif
      a.insert (i, Qxs [i]);
    }

  OsiRowCut cut;
  cut.setRow (a);

  delete [] Qxs;

  if (lambda == dCoeffLo_) {
     cut.setUb (a0);

#ifdef DEBUG
     if (optimum && (lhs - a0 > COUENNE_EPS)) {
       printf ("cut violates optimal solution: %g > %g\n", lhs, a0);
       cut.print ();
     }
     if (lhsc < a0 + COUENNE_EPS){
       printf ("cut (+) is not cutting: ");
       cut.print ();
     }
#endif
     //     cut.setLb(-COUENNE_INFINITY);
  }
  else {
    cut.setLb (a0);
#ifdef DEBUG
    if (optimum && (lhs - a0 < -COUENNE_EPS)) {
       printf ("cut violates optimal solution: %g < %g\n", lhs, a0);
       cut.print ();
    }
    if (lhsc > a0 - COUENNE_EPS){
       printf ("cut (-) is not cutting: ");
       cut.print ();
     }
#endif
    //    cut.setUb(COUENNE_INFINITY);
  }

  cs.insert (cut);
}

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


void exprQuad::quadCuts (expression *w, OsiCuts &cs, const CouenneCutGenerator *cg){

#ifdef DEBUG
  std::cout<<"Expression has "<< lcoeff_.size () <<" linear terms and "
           << nqterms_ <<" quadratic terms." << std::endl;

  printf ("Q:");
  for (sparseQ::iterator row = matrix_.begin (); row != matrix_.end (); ++row) {
    int xind = row -> first -> Index ();
    for (sparseQcol::iterator col = row -> second.begin (); col != row -> second.end (); ++col)
      printf (" <%d,%d,%g>",  xind, col -> first -> Index (), col -> second);
  }

  printf ("\nb:");
  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el)
    //for (int i=0; i < nlterms_; i++)
    printf (" <%d,%g>",  el -> first -> Index (), el -> second);//index_ [i], coeff_ [i]);

  if (c0_) 
    printf ("; <c0 = %g>", c0_);

  printf ("\nBounds: var           val        lb        ub        eigval   scaled\n");

  int index = 0;

  for (std::map <exprVar *, std::pair <CouNumber, CouNumber> >::iterator i = bounds_.begin ();
       i != bounds_.end (); ++i, index++) {

    printf ("%3d:\t", index);
    i -> first -> print (); printf ("\t");
    printf (" %8g [%8g, %8g]",
	    (*(i -> first)) (), i -> second.first, i -> second.second);

    CouNumber 
      lb = cg -> Problem () -> Lb (i -> first -> Index ()),
      ub = cg -> Problem () -> Ub (i -> first -> Index ());

    if ((eigen_.size () > 0) &&
	(fabs (ub-lb) > COUENNE_EPS))
      printf (" --> %8g %8g", 
	      eigen_.begin () -> first,
	      eigen_.begin () -> first / (ub-lb));

    printf ("\n");
  }
#endif

  // Get on which side constraint is violated to get the good lambda

  CouNumber
    varVal  = (*w)    (), 
    exprVal = (*this) (),
    lambda  =
    (eigen_.size () == 0) ? 0. :
    (varVal < exprVal) ?
      CoinMin (0., eigen_.begin  () -> first) : // Use under-estimator
      CoinMax (0., eigen_.rbegin () -> first),  // Use  over-estimator
    convVal = 0.;

  const CouenneProblem& problem = *(cg -> Problem ());
  const int numcols = problem.nVars ();

  const double
    *colsol = problem.X  (), // current solution
    *lower  = problem.Lb (), //         lower bound
    *upper  = problem.Ub (); //         upper

  // compute lower or upper convexification and check if it contains
  // the current point

  if (fabs (lambda) > COUENNE_EPS) {

    convVal = exprVal;

    // there is a convexification, check if out of current point

    for (std::map <exprVar *, std::pair <CouNumber, CouNumber> >::iterator i = bounds_.begin ();
	 i != bounds_.end (); ++i) {

      int ind = i -> first -> Index ();

      CouNumber
	xi = colsol [ind],
	lb = lower  [ind],
	ub = upper  [ind],
	delta = ub-lb;

      if (fabs (delta) > COUENNE_EPS)
	convVal += lambda * (xi-lb) * (ub-xi) / (delta * delta);
    }

    if (varVal < exprVal) {if (convVal < varVal) return;}
    else                  {if (convVal > varVal) return;}
  }

#ifdef DEBUG
  std::cout << "Point to cut: "; 
  for (int i = 0 ; i < numcols ; i++) std::cout << colsol [i] << ", ";
  printf (" (w,f(x),c) = (%g,%g,%g) -- lambda = %g\n", (*w) (), exprVal, convVal, lambda);
#endif

  // Initialize by copying $a$ into a dense vector and computing Q x^*
  double 
    *Qxs = new double [numcols], // sparse coefficient vector, $Qx^*$
     a0  = -c0_;                 // constant term

  CoinFillN (Qxs, numcols, 0.);

  // Compute 2 * Q x^*.
  for (sparseQ::iterator row = matrix_.begin (); row != matrix_.end (); ++row) {

    int qi = row -> first -> Index ();

    for (sparseQcol::iterator col = row -> second.begin (); col != row -> second.end (); ++col) {

      int qj = col -> first -> Index ();

      CouNumber 
	qc = col -> second,
	xi = colsol [qi],
	xj = colsol [qj];

      if (qi != qj) {
	Qxs [qi] += qc * xj; // contribution of element $q_{ij}$ to (Qx)_i
	Qxs [qj] += qc * xi; //                         $q_{ij}$    (Qx)_j
	a0 += 2 * qc * xi * xj;
      }
      else {
	/*
	if (fabs (lambda) > COUENNE_EPS) {

	  CouNumber
	    lb = lower  [qi],
	    ub = upper  [qi],
	    delta = ub-lb;

	  if (fabs (delta) > COUENNE_EPS)
	    qc -= lambda / (delta*delta);
	}
	*/
	// elements on the diagonal are not halved upon reading
	a0 += qc * xi * xi;
	Qxs [qi] += 2 * qc * xi; 
      }
    }
  }

#ifdef DEBUG
  printf ("2Qx = ("); for (int i = 0; i < numcols; i++) printf ("%g ", Qxs [i]); printf (")\n");
#endif

  // Add a to it.
  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el)
    Qxs [el -> first -> Index ()] += el -> second; //coeff_ [i];

  // multiply Qx^* by x^*^T again and store the result for the lower
  // bound into constant term

  /*
  for (int i=0; i < numcols; i++){
    a0 -= 0.5 * Qxs [i] * colsol [i];
    //    Qxs [i] *= 2;
  }
  */

  // And myself
  Qxs [w -> Index ()] -= 1;

#ifdef DEBUG
  printf ("2Qx = ("); for(int i = 0; i < numcols; i++) printf ("%g ", Qxs [i]); printf (")[%g]\n",a0);
#endif

  //a0 -= exprVal;

  if (fabs (lambda) > COUENNE_EPS) // Now the part which depends on lambda, if there is one

    for (std::map <exprVar *, std::pair <CouNumber, CouNumber> >::iterator i = bounds_.begin ();
	 i != bounds_.end (); ++i) {

      int ind = i -> first -> Index ();

      CouNumber 
	xi    = colsol [ind],
	lb    = lower [ind],
	ub    = upper [ind],
	delta = ub-lb;

      if (fabs (delta) > COUENNE_EPS) {

	CouNumber normlambda = lambda / (delta*delta),
	  coeff = normlambda * (lb + ub - 2. * xi);

	a0 += normlambda * (lb*ub - xi*xi);

	//a0 += coeff * xi - normlambda * (xi - lb) * (ub - xi);
	//a0 += normlambda * lb * ub;
	Qxs [ind] += coeff;
	//Qxs [ind] += normlambda * (lb + ub);
      }// else coeff = 0.;

      //      a0 += lambda [k] * lower  [ind] * upper  [ind];
      //      a0 -= lambda [k] * colsol [ind] * colsol [ind];

      //Qxs [ind] -= lambda [k] * (colsol [ind]) * 2;
    }

  // Count the number of non-zeroes
  int nnz = 0;
  for (int i=0; i < numcols ; i++)
    if (fabs (Qxs [i]) > COUENNE_EPS)
      nnz++;

#ifdef DEBUG
  printf ("2Qx = (");for(int i=0;i<numcols;i++)printf("%g ",Qxs[i]);printf (")[%g], %d nz\n",a0,nnz);
#endif

  // Pack the vector into a CoinPackedVector and generate the cut.
  CoinPackedVector a (false);
  a.reserve (nnz);

#ifdef DEBUG
  CouNumber lhs = 0, lhsc = 0,
    *optimum = cg -> Problem () -> bestSol (),
    *current = cg -> Problem () -> X ();
#endif

  for (int i=0; i < numcols; i++)

    if (fabs (Qxs [i]) > 1.0e-21) { // why 1.0e-21? Look at CoinPackedMatrix.cpp:2188

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
    
  if (varVal < exprVal) { //(lambda == dCoeffLo_) {

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

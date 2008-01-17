/*
 * Name:    minMaxDelta.cpp
 * Author:  Pietro Belotti
 * Purpose: general function for computing best branching point based
 *          on min max height of resulting convexifications (dychotomic 
 *          search)
 *
 * (C) Carnegie-Mellon University, 2007.
 * This file is licensed under the Common Public License (CPL)
 */

//#define TEST_BRPTS
/*#ifdef  TEST_BRPTS

#include <stdio.h>
#include <stdlib.h>

#include "CouennePrecisions.hpp"
typedef double CouNumber;
typedef CouNumber (*unary_function) (CouNumber);
CouNumber midInterval (CouNumber, CouNumber, CouNumber);
CouNumber safe_pow (register CouNumber base, 
		    register CouNumber exponent);

///
class funtriplet {

public:

  /// Basic constructor
  funtriplet () {}

  /// Destructor
  virtual ~funtriplet () {}

  virtual CouNumber F     (CouNumber x) = 0; //< main funtion
  virtual CouNumber Fp    (CouNumber x) = 0; //< first-order derivative of main funtion
  virtual CouNumber Fpp   (CouNumber x) = 0; //< second-order derivative of main funtion
  virtual CouNumber FpInv (CouNumber x) = 0; //< inverse of the first-order derivative
};


///
class simpletriplet: public funtriplet {

protected:

  unary_function f_;   //< the function 
  unary_function fp_;  //< the first-order derivative
  unary_function fpp_; //< the second-order derivative 
  unary_function fpI_; //< the inverse of the first-order derivative 

public:

  /// Basic constructor
  simpletriplet (unary_function f   = NULL, 
		 unary_function fp  = NULL, 
		 unary_function fpp = NULL,
		 unary_function fpI = NULL):
    f_   (f),
    fp_  (fp),
    fpp_ (fpp),
    fpI_ (fpI) {}

  /// Destructor
  virtual ~simpletriplet () {}

  virtual CouNumber F     (CouNumber x) {return f_   (x);} //< main funtion
  virtual CouNumber Fp    (CouNumber x) {return fp_  (x);} //< first-order derivative
  virtual CouNumber Fpp   (CouNumber x) {return fpp_ (x);} //< second-order derivative
  virtual CouNumber FpInv (CouNumber x) {return fpI_ (x);} //< inverse of first-order derivative
};


///
class powertriplet: public funtriplet {

protected:

  CouNumber exponent_; //< defines the power function triplet

public:

  /// Basic constructor
  powertriplet (CouNumber exponent):
    exponent_ (exponent) {}

  /// Destructor
  virtual ~powertriplet () {}

  virtual CouNumber F   (CouNumber x) 
  {return safe_pow (x, exponent_);}                                   //< main funtion

  virtual CouNumber Fp  (CouNumber x) 
  {return exponent_ * safe_pow (x, exponent_ - 1);}                   //< first-order derivative 

  virtual CouNumber Fpp (CouNumber x) 
  {return exponent_ * (exponent_ - 1) * safe_pow (x, exponent_ - 2);} //< second-order derivative 

  virtual CouNumber FpInv (CouNumber x) 
  {return safe_pow (x / exponent_, 1 / (exponent_ - 1));} //< inverse of first derivative
};


///
class kpowertriplet: public powertriplet {

protected:

  CouNumber mult_; //< pre-multiplier

public:

  /// Basic constructor
  kpowertriplet (CouNumber exponent, CouNumber k):
    powertriplet (exponent),
    mult_ (k) {}

  /// Destructor
  virtual ~kpowertriplet () {}

  virtual CouNumber F   (CouNumber x)  //< main funtion
  {return mult_ * safe_pow (x, exponent_);}

  virtual CouNumber Fp  (CouNumber x)  //< first-order derivative 
  {return mult_ * exponent_ * safe_pow (x, exponent_ - 1);}

  virtual CouNumber Fpp (CouNumber x)  //< second-order derivative 
  {return mult_ * exponent_ * (exponent_ - 1) * safe_pow (x, exponent_ - 2);}

  virtual CouNumber FpInv (CouNumber x) 
  {return safe_pow (x / (mult_ * exponent_), 1 / (exponent_ - 1));} //< inverse of first derivative
};



#else */
#include "CouenneObject.hpp"
#include "funtriplets.hpp"
//#endif


const int maxIter = 20;


void choosePoint (CouNumber x0, CouNumber y0, 
		  CouNumber l, CouNumber u, 
		  CouNumber &brpt)   {

  double nl, nu;

  if (x0 < brpt) {
    nl = l;
    nu = brpt;
  } else {
    nl = brpt;
    nu = u;
  }

  double
    lnl = log (nl),
    lnu = log (nu),
    asq = lnu-lnl,
    bsq = u-l,
    distslope = fabs ((lnu-lnl) * (x0-nl) - (nu-nl) * (y0-lnl)) /
    sqrt (asq*asq + bsq*bsq),
    csq = x0 - exp (y0),
    dsq = y0 - log (x0),
    distcurve = 0.5 * sqrt (csq*csq + dsq*dsq);

  if (distcurve < distslope)
    brpt = midInterval (x0, l, u); 
}


///
CouNumber curvDistance (funtriplet *ft, CouNumber lb, CouNumber ub) {

  // Consider the function f(x) between lb and ub. The slope of the
  // convexification on the concave side, y = alpha x + alpha0, is:

  CouNumber alpha = (ft -> F (ub) - ft -> F (lb)) / (ub - lb);

  // and the constant term, alpha0, is

  CouNumber alpha0 = (ub * ft -> F (lb) - lb * ft -> F (ub)) / (ub - lb);

  // The point at which f(.) has derivative equal to the slope is the
  // point of maximum height w.r.t the slope. The point z where
  // maximum of f(z) - (ax+b), where (ax+b) is the convexification
  // line (a=slope), is such that 
  //
  // f'(z) - alpha = 0   ==> z = (f')^{-1} (alpha)

  CouNumber z = ft -> FpInv (alpha);

  // The real height is computed as [f(z) - (alpha * z + alpha0)]
  // divided by the norm sqrt (alpha^2 + 1)

  return ((ft -> F (z) - (alpha * z + alpha0)) / sqrt (alpha * alpha + 1));
}


///
CouNumber minMaxDelta (funtriplet *ft, 
		       CouNumber x, CouNumber y, 
		       CouNumber lb, CouNumber ub) {

  CouNumber 
    lbm = lb,                // extremes of the interval where to look 
    ubm = ub,     
    b   = 0.5 * (lbm + ubm); // starting point

  for (int iter = 0; iter < maxIter; iter++) {

    CouNumber distL = curvDistance (ft, lb,  b),  // max height at left
              distR = curvDistance (ft,  b, ub),  // max height at right
              delta = fabs (distL) - fabs (distR);

    //    fprintf (stderr, "%4d %10g %10g %10g %10g %10g %10g\n", 
    //	     iter, lbm, ubm, b, distL, distR, delta);

    if (fabs (delta) < COUENNE_EPS) 
      break;

    CouNumber oldb = b;

    // choose a smarter b based on an estimate of the derivative of
    // the distance function at the current point, knowing it's null
    // at the extremes

    /*
    if ((iter < 0) ||
	(distL * distR < 0)            || // weird, they have opposite sign
	(fabs (b - lbm) < COUENNE_EPS) || // too close to lower bound
	(fabs (b - ubm) < COUENNE_EPS))   // too close to upper bound
    */

      b = 0.5 * (lbm + ubm);

    /*
    else 

      b = (distR * ubm / (ubm - b) + 
	   distL * lbm / (b - lbm)) /
	  (distL / (b - lbm) + 
	   distR / (ubm - b));
    */

    if (delta > 0) ubm = oldb; // right max height is smaller, move left
    else           lbm = oldb; // and viceversa
  }

  return midInterval (b, lb, ub);
}


///
CouNumber maxHeight (funtriplet *ft, 
		     CouNumber x, CouNumber y, 
		     CouNumber lb, CouNumber ub) {
  /* fprintf (stderr,"slope is (%g - %g) / (%g - %g) = %g / %g = %g ----> inverse is %g\n", 
	  ft -> F (ub), 
	  ft -> F (lb), 
	  ub, lb,
	  ft -> F (ub) - ft -> F (lb),
	  (ub - lb),
	  (ft -> F (ub) - ft -> F (lb)) / (ub - lb),
	  ft -> FpInv ((ft -> F (ub) - ft -> F (lb)) / (ub - lb)));*/
  return (ft -> FpInv ((ft -> F (ub) - ft -> F (lb)) / (ub - lb)));
}



// testing stuff ////////////////////////////////////////////////////////////////////////

/*#ifdef TEST_BRPTS

inline CouNumber safe_pow (register CouNumber base, 
			   register CouNumber exponent) {

  if (base < 0) {

    register int rndexp;

    if (((fabs (exponent - (rndexp = COUENNE_round (exponent))) < COUENNE_EPS) ||
	 ((fabs (exponent) > COUENNE_EPS) && 
	  (fabs (1. / exponent - (rndexp = COUENNE_round (1. / exponent))) < COUENNE_EPS)))
	&& (rndexp % 2))
      return (- pow (- base, exponent));
  }

  if (fabs (base) >= COUENNE_INFINITY) {

    if (base <= -COUENNE_INFINITY) {

      register int intk = COUENNE_round (exponent);

      if ((fabs (exponent - intk) < COUENNE_EPS) && (intk % 2))
	return (exponent < 0) ? 0 : -COUENNE_INFINITY;
    }
    else return (exponent < 0) ? 0 : COUENNE_INFINITY;
  }

  return (pow (base, exponent));
}


CouNumber midInterval (CouNumber curr, CouNumber l, CouNumber u) {

  CouNumber x = curr;

  if      (x<l) x = l;
  else if (x>u) x = u;

#define COUENNE_NEAR_BOUND 1e-4

  CouNumber alpha = 0.25, retval;//CouenneBranchingObject::Alpha (), retval;
 
  if ((x-l < COUENNE_NEAR_BOUND) ||
      (u-x < COUENNE_NEAR_BOUND))
    if      (u <   COUENNE_INFINITY)
      if    (l > - COUENNE_INFINITY) retval = alpha * x + (1. - alpha) * (l + u) / 2.;
      else                           retval = -1 + (u<0) ? u*2 : u/2;
    else if (l > - COUENNE_INFINITY) retval = +1 + (l>0) ? l*2 : l/2;
    else                             retval = 0;
  else retval = x;//alpha * x + (1. - alpha) * (l + u) / 2.;

  return retval;
}

CouNumber inv (CouNumber arg) {return 1./arg;}
CouNumber negInvSqr (CouNumber arg) {return -1./(arg*arg);}


int main (int argc, char **argv) {

  CouNumber 
    lb = atof (argv [1]),
    ub = atof (argv [2]), b;

  //  simpletriplet ft (exp, exp, exp, log); // exponential 
  //  simpletriplet ft (log, inv, negInvSqr, inv); // logarithm
  powertriplet ft (4);

#define NPTS 500

  for (int i=0; i <= NPTS; i++) {
    double iter = lb + (ub-lb) / (double) NPTS * (double) i;
    printf ("%g %g\n", iter, ft. F (iter));
  }

  printf ("%g %g\n", lb, ft. F (lb));

  b = minMaxDelta (&ft, b, 0, lb, ub);
  printf ("%g %g\n", b, ft. F (b));
  printf ("%g %g\n", ub, ft. F (ub));

  fprintf (stderr, "mmd = %g\n", b);

  b = maxHeight (&ft, 1, 0, lb, ub);
  printf ("%g %g\n", b, ft. F (b));
  //printf ("%g %g\n", b, ft. F (b) - 1);
  //printf ("%g %g\n", b, ft. F (b));
  printf ("%g %g\n", lb, ft. F (lb));

  fprintf (stderr, "maxHeight = %g\n", b);
}
#endif*/

/*
 * Name:    impliedBounds-sum.cpp
 * Author:  Pietro Belotti
 * Purpose: inferring bounds on monomials in a sum
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>

#include "CoinHelperFunctions.hpp"
#include "exprSum.hpp"

#define ALMOST_INF (1e-5 * COUENNE_INFINITY)

int exprSum::impliedBoundSum (CouNumber wl, 
			      CouNumber wu, 
			      std::vector <CouNumber> &xl,
			      std::vector <CouNumber> &xu,
			      std::vector <std::pair <int, CouNumber> > &nl,
			      std::vector <std::pair <int, CouNumber> > &nu) {
  CouNumber 
    sumLB = 0,
    sumUB = 0;

  int 
    nImpr = 0,
    n     = xl.size (), 
    infLo = -1, 
    infUp = -1;

  // check lower bounds
  for (int i=0; i<n; i++) {
    CouNumber l = xl [i];
    if (l < -ALMOST_INF) 
      if (infLo >= 0) {infLo = -2; break;}
      else infLo = i;
    else sumLB += l;
  }

  // check upper bounds
  for (int i=0; i<n; i++) {
    CouNumber u = xu [i];
    if (u > ALMOST_INF) 
      if (infUp >= 0) {infUp = -2; break;}
      else infUp = i;
    else sumUB += u;
  }

  // if more than two unbounded quantities on both sides, bail out
  if ((infUp == -2) && 
      (infLo == -2)) 
    return 0;

  // new lower bounds ////////////////////////////////////////////////////

  if (infLo == -1) { // none of the "first" components is unbounded from below

    for (int i=0; i<n; i++) {

      CouNumber nb = wu - sumLB + xl [i];
      if (nb < xu [i]) {
	nu.push_back (std::pair <int, CouNumber> (i, nb));
	nImpr ++;
      }
    }	

  } else if (infLo >= 0) { // exactly one of them is, can improve bound on that one only
  
    CouNumber nb = wu - sumLB;
    if (nb < xu [infLo]) {
      nu.push_back (std::pair <int, CouNumber> (infLo, nb));
      nImpr ++;
    }
  }

  // new upper bounds ////////////////////////////////////////////////////
 
  if (infUp == -1) { // none of the "first" components is unbounded from below

    for (int i=0; i<n; i++) {

      CouNumber nb = wl - sumUB + xu [i];
      if (nb > xl [i]) {
	nl.push_back (std::pair <int, CouNumber> (i, nb));
	nImpr ++;
      }
    }	

  } else if (infUp >= 0) { // exactly one of them is, can improve bound on that one only
  
    CouNumber nb = wl - sumUB;
    if (nb < xl [infLo]) {
      nl.push_back (std::pair <int, CouNumber> (infUp, nb));
      nImpr ++;
    }
  }

  return nImpr;
}

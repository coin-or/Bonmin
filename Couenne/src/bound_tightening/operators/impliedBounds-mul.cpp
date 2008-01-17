/*
 * Name:    impliedBounds-mul.cpp
 * Author:  Pietro Belotti
 * Purpose: inferring bounds on factors of a product
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>

#include "CoinHelperFunctions.hpp"
#include "exprMul.hpp"

int exprMul::impliedBoundMul (CouNumber wl, 
			      CouNumber wu, 
			      std::vector <CouNumber> &xlv,
			      std::vector <CouNumber> &xuv,
			      std::vector <std::pair <int, CouNumber> > &nl,
			      std::vector <std::pair <int, CouNumber> > &nu) {
  int nImpr = 0;

  // case 1: vectors of order 2

  if (xlv.size () == 2) {

    CouNumber 
      xl = xlv [0],  xu = xuv [0],
      yl = xlv [1],  yu = xuv [1];

    if (wl >= 0.) {

      // point B in central infeasible area

      if (xu * yu < wl) {
	if (xu * yl < wl) {nu.push_back (std::pair <int, CouNumber> (0, wl / yl)); nImpr++;}
	if (xl * yu < wl) {nu.push_back (std::pair <int, CouNumber> (1, wl / xl)); nImpr++;}
	//resxU = (*xu * *yl < wl) && updateBound (+1, xu, wl / *yl);
	//resyU = (*xl * *yu < wl) && updateBound (+1, yu, wl / *xl);
      }

      // point C in central infeasible area

      if (xl * yl < wl) {
	if (xl * yu < wl) {nl.push_back (std::pair <int, CouNumber> (0, wl / yu)); nImpr++;}
	if (xu * yl < wl) {nl.push_back (std::pair <int, CouNumber> (1, wl / xu)); nImpr++;}
	//resxL = (xl * yu < wl) && updateBound (-1, xl, wl / yu);
	//resyL = (xu * yl < wl) && updateBound (-1, yl, wl / xu);
      }
    } else if (wl > -COUENNE_INFINITY) {

      // the infeasible set is a hyperbola with two branches

      // upper left
      if ((xl*yl < wl) && (yl>0.)) {nl.push_back (std::pair <int, CouNumber> (0, wl / yl)); nImpr++;}
      if ((xu*yu < wl) && (yu>0.)) {nu.push_back (std::pair <int, CouNumber> (1, wl / xu)); nImpr++;}
      //resxL = (xl * yl < wl) && (yl > 0.) && updateBound (-1, xl, wl / yl); // point C
      //resyU = (xu * yu < wl) && (yu > 0.) && updateBound (+1, yu, wl / xu); // point B

      // lower right
      if ((xu*yu < wl) && (yu<0.)) {nu.push_back (std::pair <int, CouNumber> (0, wl / yu)); nImpr++;}
      if ((xl*yl < wl) && (yl<0.)) {nl.push_back (std::pair <int, CouNumber> (1, wl / xl)); nImpr++;}
      //resxU = (xu * yu < wl) && (yu < 0.) && updateBound (+1, xu, wl / yu); // point B
      //resyL = (xl * yl < wl) && (yl < 0.) && updateBound (-1, yl, wl / xl); // point C
    }


    // w's upper bound ///////////////////////////////////////////

    if (wu >= 0.) {

      if (wu < COUENNE_INFINITY) {
	// the infeasible set is a hyperbola with two branches

	// upper right
	if ((xu*yl > wu) && (yl>0.)) {nu.push_back (std::pair <int, CouNumber> (0, wu/yl)); nImpr++;}
	if ((xl*yu > wu) && (yu>0.)) {nu.push_back (std::pair <int, CouNumber> (1, wu/xl)); nImpr++;}
	//resxU = (xu * yl > wu) && (yl > 0.) && updateBound (+1, xu, wu / yl) || resxU; // point D
	//resyU = (xl * yu > wu) && (yu > 0.) && updateBound (+1, yu, wu / xl) || resyU; // point A

	// lower left
	if ((xl*yu > wu) && (yu<0.)) {nl.push_back (std::pair <int, CouNumber> (0, wu/yu)); nImpr++;}
	if ((xu*yl > wu) && (yl<0.)) {nl.push_back (std::pair <int, CouNumber> (1, wu/xu)); nImpr++;}
	//resxL = (xl * yu > wu) && (yu < 0.) && updateBound (-1, xl, wu / yu) || resxL; // point A
	//resyL = (xu * yl > wu) && (yl < 0.) && updateBound (-1, yl, wu / xu) || resyL; // point D
      }

    } else {

      // point D in central infeasible area

      if (xu * yl > wu) {
	if (xu * yu > wu) {nu.push_back (std::pair <int, CouNumber> (0, wu / yu)); nImpr++;}
	if (xl * yl > wu) {nl.push_back (std::pair <int, CouNumber> (1, wu / xl)); nImpr++;}
	//resxU = (xu * yu > wu) && updateBound (+1, xu, wu / yu) || resxU;
	//resyL = (xl * yl > wu) && updateBound (-1, yl, wu / xl) || resyL;
      }

      // point A in central infeasible area

      if (xl * yu > wu) {
	if (xl * yl > wu) {nl.push_back (std::pair <int, CouNumber> (0, wu / yl)); nImpr++;}
	if (xu * yu > wu) {nu.push_back (std::pair <int, CouNumber> (1, wu / xu)); nImpr++;}
	//resxL = (xl * yl > wu) && updateBound (-1, xl, wu / yl) || resxL;
	//resyU = (xu * yu > wu) && updateBound (+1, yu, wu / xu) || resyU;
      }
    }
  }

  return nImpr;
}

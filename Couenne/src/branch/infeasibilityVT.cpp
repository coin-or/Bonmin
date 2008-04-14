/*
 * Name:    infeasibilityVT.cpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: Compute violation transfer of a variable
 *
 * (C) Carnegie-Mellon University, 2008.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinHelperFunctions.hpp"

#include "CouenneProblem.hpp"
#include "CouenneVTObject.hpp"

const CouNumber weiMin = 0.3;
const CouNumber weiMax = 1.3;
const CouNumber weiSum = 1.1; // at least 1 so top level aux are avoided
const CouNumber weiAvg = 0.0;


/// compute infeasibility of this variable, |w - f(x)| (where w is
/// the auxiliary variable defined as w = f(x))
double CouenneVTObject::infeasibility (const OsiBranchingInformation *info, int &way) const {

  assert (reference_);
  int index = reference_ -> Index ();

  assert (index >= 0);
  const std::set <int> &dependence = problem_ -> Dependence () [index];

  assert (reference_ -> Image ());
  CouNumber retval = fabs (((*reference_ -> Image ())) () - (*reference_) ());

  if (dependence.size () == 0) { // this is a top level auxiliary,
				 // nowhere an independent

    // that means, for VT, that letting w vary and keeping everything
    // else constant we return the difference between current value
    // and value of function at this point

    if (retval < CoinMin (COUENNE_EPS, feas_tolerance_)) 
      retval = 0;

  } else { // this appears as an independent in all auxs of the
	   // "dependence" list

    CouNumber 
      xcurr = info -> solution_ [index],
      lFeas = xcurr, 
      rFeas = xcurr;

    for (std::set <int>::const_iterator i = dependence.begin ();
	 i != dependence.end (); ++i) {

      CouNumber 
	left  = xcurr, 
	right = xcurr;

      const CouenneObject &obj = problem_ -> Objects () [*i];
      if (obj. Reference () && (retval >= CoinMin (COUENNE_EPS, feas_tolerance_))) {

	obj.Reference () -> extendInterval (info, reference_, left, right);

	if (left  < lFeas)  lFeas  = left;
	if (right > lRight) lRight = right;
      }
    }

    if (lFeas < info -> lower_ [index]) lFeas = info -> lower_ [index];
    if (rFeas > info -> upper_ [index]) rFeas = info -> upper_ [index];

    retval = rFeas - lFeas;

    const threshold = .2;

    if (retval > COUENNE_EPS) {

      CouNumber
	leanLeft  = (xcurr - lFeas) / retval;

      if      (leanLeft <     threshold) way = 0;
      else if (leanLeft > 1 - threshold) way = 1;
    }
  }

  // done computing delta. Now transfer violation on LP relaxation
  //
  // info -> pi_                         duals
  // info -> solution_                   solution
  // info -> lower_                      lower bound
  // info -> upper_                      upper bound
  // info -> solver () -> getNumRows()   # rows
  // info -> numberColumns_              # cols

  CouNumber vt_delta = 0;

  if (reference_ -> Index () == problem_ -> Obj (0) -> Body () -> Index ())
    vt_delta += 1;

  int 
    nRows = info -> solver () -> getNumRows (),
    nCols = info -> numberColumns_;

  for (int i=0; i<nRows; i++)
    vt_delta += 0;

  retval = alpha * vt_delta + beta * retval;

  if ((retval > CoinMin (COUENNE_EPS, feas_tolerance_)) &&
      (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING))) {

    printf ("vt-delta is %-10g [", retval); 

    reference_ -> print (); 
    if (dependence.size () == 0) { // if no list, print image
      printf (" := ");
      reference_ -> Image () -> print ();
    }
    printf ("]\n");
  }

  // add term to stop branching on very tiny intervals

  return ((retval < CoinMin (COUENNE_EPS, feas_tolerance_)) ? 
	  0. : retval);
}

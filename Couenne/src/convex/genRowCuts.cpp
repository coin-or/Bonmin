/*
 * Name:    genRowCuts.cpp
 * Author:  Pietro Belotti
 * Purpose: generate Row Cuts for current convexification
 *
 * (C) Pietro Belotti, all rights reserved. 
 * This file is licensed under the Common Public License.
 */

#include <CglCutGenerator.hpp>
#include <CouenneCutGenerator.h>
#include <CouenneProblem.h>


/// generate OsiRowCuts for current convexification
void CouenneCutGenerator::genRowCuts (const OsiSolverInterface &si,
				      OsiCuts &cs,
				      int nchanged, 
				      int *changed,
				      t_chg_bounds *chg) const {

  // For each auxiliary variable, create convexification cut (or set
  // of cuts) and add it to cs

  if (firstcall_)
    for (int i=0, j = problem_ -> nAuxs (); j--; i++) {
      if (problem_ -> Aux (i) -> Multiplicity () > 0)
	problem_ -> Aux (i) -> generateCuts (si, cs, this, chg);
    }
  else { // chg_bds contains the indices of the variables whose bounds
	 // have changes (a -1 follows the last element)

    for (int i = 0, j = problem_ -> nAuxs (); j--; i++) {

      // TODO: check if list contains all and only aux's to cut

      expression * image = problem_ -> Aux (i) -> Image ();

      // cut only if:

      if (   (image -> Linearity () > LINEAR)         // 1) expression is non linear
	  && (image -> dependsOn (changed, nchanged)) // 2) it depends on changed variables
	     || !changed)                             // or, if no changed variable is passed

      if (problem_ -> Aux (i) -> Multiplicity () > 0)
	problem_ -> Aux (i) -> generateCuts (si, cs, this, chg);
    }
  }
}

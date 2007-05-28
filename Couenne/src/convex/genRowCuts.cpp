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
				      int *changed) const {

  // For each auxiliary variable, create cut (or set of cuts) violated
  // by current point and add it to cs

  if (firstcall_)
    for (int i=0, j = problem_ -> nAuxs (); j--; i++) {
      if (problem_ -> Aux (i) -> Multiplicity () > 0)
	problem_ -> Aux (i) -> generateCuts (si, cs, this);
    }
  else { // chg_bds contains the indices of the variables whose bounds
	 // have changes (a -1 follows the last element)

    for (int i=0, j = problem_ -> nAuxs (); j--; i++) {

      // TODO: check if list contains all and only aux's to cut
      /*expression * image = problem_ -> Aux (i) -> Image ();
      if (   (image -> Linearity () > LINEAR)          // if linear, no need to cut twice

	  && (image -> dependsOn (changed, nchanged)) // if expression does not depend 
                                                      // on changed variables, do not cut 
          || !changed)*/                              // or if no changed variable is passed
      if (problem_ -> Aux (i) -> Multiplicity () > 0)
	problem_ -> Aux (i) -> generateCuts (si, cs, this);
    }
  }
}

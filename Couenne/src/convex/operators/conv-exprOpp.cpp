/*
 * Name:    conv-exprOpp.cpp
 * Author:  Pietro Belotti
 * Purpose: methods to convexify opposite of expressions
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprOpp.hpp>
#include <exprConst.hpp>

#include <CouenneProblem.hpp>
#include <CouenneCutGenerator.hpp>


// generate equality between *this and *w
void exprOpp::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg,
			    t_chg_bounds *chg, 
			    int wind, CouNumber lb, CouNumber ub) {

  // If wind = -1 then this is the normal procedure (see below,
  // "easy..."). Otherwise, there is a constraint of the form -x <= b
  // that was translated into auxiliary w = -x and w <= b. wind is w's
  // index (which we have to single out)

  if (wind >= 0) {

    int xind = argument_ -> Index ();

    if (xind < 0) {
      printf ("#### invalid index for exprOpp::gencuts()\n");
      return;
    }

    OsiColCut *cut = new OsiColCut;

    CouNumber 
      &xlb = cg -> Problem () -> Lb (xind),
      &xub = cg -> Problem () -> Ub (xind);

    if (-ub > xlb) xlb = -ub;
    if (-lb < xub) xub = -lb;

    cut -> setLbs (1, &xind, &xlb);
    cut -> setUbs (1, &xind, &xub);

    cs.insert (cut);

    delete (cut);
  }
  else // easy... 
    if (cg -> isFirst ())
      cg -> createCut (cs, 0., 0, w->Index (), 1., argument_->Index (), 1.);
}

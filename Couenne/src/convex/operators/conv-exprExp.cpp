/*
 * Name:    conv-exprExp.C
 * Author:  Pietro Belotti
 * Purpose: convexification and bounding methods for the exponential operator
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprExp.h>
#include <exprConst.h>
#include <exprAux.h>

#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>

#define MULT_FACTOR 10


// generate convexification cut for constraint w = this

void exprExp::generateCuts (exprAux *aux, const OsiSolverInterface &si, 
			    OsiCuts &cs,  const CouenneCutGenerator *cg) {

  expression *le, *ue;

  argument_ -> getBounds (le, ue);

  CouNumber w = (*aux) (),
            x = (*argument_) (),
            l = (*le) (),
            u = (*ue) ();

  bool check = cg -> isFirst () || !(cg -> addViolated ());

  int w_ind = aux       -> Index (),
      x_ind = argument_ -> Index ();

  OsiRowCut *cut;

  // upper segment

  CouNumber expl = exp (l);
  CouNumber oppslope = (expl - exp (u)) / (u - l);

  if ((u < COUENNE_INFINITY - 1) 
      && (check || 
	  (w + oppslope*x > expl + oppslope*l)) &&
      ((cut = cg -> createCut (expl + oppslope*l, -1, 
			       w_ind, CouNumber (1.), 
			       x_ind, oppslope)))) {

    //    printf ("Exp upper: "); cut -> print ();
    cs.insert (cut);
  }

  // lower convexification: start with trivial envelope w >= 0

  if (check || (w < - COUENNE_EPS)) {

    if ((cut = cg -> createCut (CouNumber (0.), +1, w_ind, CouNumber (1.)))) {

      //      printf ("Exp trivial: "); cut -> print ();
      cs.insert (cut);
    }
  }

  // add tangent points: first choose sampling points

  int ns = cg -> nSamples ();

  // fix bounds to get finite coefficients

  if (l < - COUENNE_INFINITY + 1) {
    if (u > COUENNE_INFINITY - 1) {
      
      l = - log (MULT_FACTOR) * ns/2;
      u =   log (MULT_FACTOR) * ns/2;
    } else 
      l = u - log (2) * ns;

  } else if (u > COUENNE_INFINITY - 1) 
    u = l + log (2) * ns;

  // approximate the exponential function from below
  cg -> addEnvelope (cs, +1, exp, exp, w_ind, x_ind, x, l, u);

  delete le;
  delete ue;
}

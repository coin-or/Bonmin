/*
 * Name:    conv-exprCos.C
 * Author:  Pietro Belotti
 * Purpose: convexification methods for the cosine class
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>

#include <exprCos.h>
#include <exprClone.h>
#include <exprAux.h>
#include <exprSum.h>
#include <exprSub.h>

#include <CouenneCutGenerator.h>

// construct under-estimator for cosine expression w = cos x

int exprCos::lowerLinearHull (exprAux *w, int *&nterms, expression ***&coeff, 
			      int **&indices, expression **&rhs, enum con_sign *&sign) {

  // For cosine, it suffices -- at least for now -- to add three
  // halfspaces:
  //
  // 1)  w     >= -1
  // 2)  w + x >= cos lb + lb 
  // 3)  w - x >= cos ub - ub
  //
  // the first is tangent if the bound interval includes points where
  // cosine is minimum, while 2) and 3) are simple, but not tight
  // bounds (numerical solution of an equation is needed for that)

  nterms = new int [3];
  nterms [0] = 1;
  nterms [1] = 2;
  nterms [2] = 2;

  allocateCon (3, nterms, coeff, indices, rhs, sign);

  // constraint 1)

  coeff   [0] [0] = new exprConst (1);   indices [0] [0] = w -> Index ();
  rhs     [0]     = new exprConst (-1);
  sign    [0]     = COUENNE_GE;

  expression *lb, *ub;
  argument_ -> getBounds (lb, ub);

  // constraint 2)

  coeff   [1] [0] = new exprConst (1);   indices [1] [0] = w -> Index ();
  coeff   [1] [1] = new exprConst (1);   indices [1] [1] = argument_ -> Index ();
  rhs     [1]     = new exprSum (new exprCos (lb), new exprClone (lb));
  sign    [1]     = COUENNE_GE;

  // constraint 3)

  coeff   [2] [0] = new exprConst (1);   indices [2] [0] = w -> Index ();
  coeff   [2] [1] = new exprConst (-1);  indices [2] [1] = argument_ -> Index ();
  rhs     [2]     = new exprSub (new exprCos (ub), new exprClone (ub));
  sign    [2]     = COUENNE_GE;

  return 3;
}


// similarly, construct over-estimator for expression within problem
// *p (p is used to add convexification constraints). It is also
// used when this function appears with a minus sign in the
// expression

int exprCos::upperLinearHull (exprAux *w, int *&nterms, expression ***&coeff, 
			      int **&indices, expression **&rhs, enum con_sign *&sign) {

  // Similarly, it suffices (at least for now) to add three
  // halfspaces:
  //
  // 1)  w     <= 1
  // 2)  w + x <= cos ub + ub 
  // 3)  w - x <= cos lb - lb
  //
  // the first is tangent if the bound interval includes points where
  // cosine is minimum, while 2) and 3) are simple, but not tight
  // bounds (numerical solution of an equation is needed for that)

  nterms = new int [3];
  nterms [0] = 1;
  nterms [1] = 2;
  nterms [2] = 2;

  allocateCon (3, nterms, coeff, indices, rhs, sign);

  // constraint 1)

  coeff   [0] [0] = new exprConst (1);   indices [0] [0] = w -> Index ();
  rhs     [0]     = new exprConst (1);
  sign    [0]     = COUENNE_LE;

  expression *lb, *ub;
  argument_ -> getBounds (lb, ub);

  // constraint 2)

  coeff   [1] [0] = new exprConst (1);   indices [1] [0] = w -> Index ();
  coeff   [1] [1] = new exprConst (1);   indices [1] [1] = argument_ -> Index ();
  rhs     [1]     = new exprSum (new exprCos (ub), new exprClone (ub));
  sign    [1]     = COUENNE_LE;                                             

  // constraint 3)

  coeff   [2] [0] = new exprConst (1);   indices [2] [0] = w -> Index ();
  coeff   [2] [1] = new exprConst (-1);  indices [2] [1] = argument_ -> Index ();
  rhs     [2]     = new exprSub (new exprCos (lb), new exprClone (lb));
  sign    [2]     = COUENNE_LE;

  return 3;
}


// generate convexification cut for constraint w = this

void exprCos::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  trigGenCuts (w, cs, cg, cos);
}

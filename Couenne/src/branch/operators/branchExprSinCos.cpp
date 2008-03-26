/*
 * Name:    branchExprSinCos.cpp
 * Author:  Pietro Belotti
 * Purpose: return branch gain and branch object for sines/cosines
 *
 * (C) Carnegie-Mellon University, 2006-07.
 * This file is licensed under the Common Public License (CPL)
 */

#include <math.h>

#include "funtriplets.hpp"
#include "exprSin.hpp"
#include "CouenneObject.hpp"
#include "CouenneBranchingObject.hpp"

static inline double oppcos  (double x) {return -cos (x);}
static inline double oppsin  (double x) {return -sin (x);}
static inline double oppasin (double x) {return asin (-x);}


/// generalized procedure for both sine and cosine
CouNumber trigSelBranch (const CouenneObject *obj, 
			 const OsiBranchingInformation *info,
			 expression *&var,
			 double * &brpts, 
			 int &way,
			 enum cou_trig type) {

  exprVar *ref = obj -> Reference ();

  var = ref -> Image () -> Argument ();

  assert (var -> Index () >= 0);
  assert (ref -> Index () >= 0);

  CouNumber l, u,
    x0 = info -> solution_ [var -> Index ()],
    y0 = info -> solution_ [ref -> Index ()];

  var -> getBounds (l,u);
				 
  simpletriplet ft ((type == COU_SINE) ? (unary_function) sin    : (unary_function) cos, 
		    (type == COU_SINE) ? (unary_function) cos    : oppsin, 
		    (type == COU_SINE) ?                  oppsin : oppcos, 
		    (type == COU_SINE) ? (unary_function) acos   : oppasin);

  brpts = (double *) realloc (brpts, sizeof (double));
  *brpts = obj -> getBrPoint (&ft, x0, l, u);

  return (y0 - ((type == COU_SINE) ? sin(x0) : cos(x0)));
}

/*
 * Name:    fake_tightening.cpp
 * Author:  Pietro Belotti
 * Purpose: fake single bounds in variables to exclude parts of the solution space 
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "BonBabInfos.hpp"

//#define DEBUG

#define MAX_ITER 10 // max # fake tightening (inner) iterations 
#define AGGR_MUL  2 // the larger,  the more conservative. Must be > 0
#define AGGR_DIV  2 // the smaller, the more conservative. Must be > 1

const CouNumber phi = 0.5 * (1. + sqrt (5.));

// core of the bound tightening procedure
bool btCore (const CouenneCutGenerator *cg,
	     const OsiSolverInterface *psi,
	     OsiCuts &cs, 
	     t_chg_bounds *chg_bds, 
	     Bonmin::BabInfo * babInfo,
	     bool serious);


// create fictitious bounds to tighten current interval
CouNumber fictBounds (char direction,
		      CouNumber  x,
		      CouNumber  lb,   
		      CouNumber  ub) {

  if   (lb < -COUENNE_INFINITY / 10) {
    if (ub >  COUENNE_INFINITY / 10) { // ]-inf,+inf[

      if (fabs (x) < COUENNE_EPS) return (direction ? AGGR_MUL : - AGGR_MUL);
      else                        return (direction ? AGGR_MUL : - AGGR_MUL) * fabs (x);

    } else { // ]-inf,u]

      if      (x < -COUENNE_EPS) return (direction ? CoinMin (0., (x+ub)/2) : AGGR_MUL * x);
      else if (x >  COUENNE_EPS) return (direction ? (x + (ub-x)/AGGR_DIV) : 0);
      else                       return (direction ? (ub/AGGR_DIV) : -AGGR_MUL);
    }
  }
  else {
    if (ub >  COUENNE_INFINITY / 10) { // [l,+inf[

      if      (x < -COUENNE_EPS) return (direction ? 0 : (x - (x-lb) / AGGR_DIV));
      else if (x >  COUENNE_EPS) return (direction ? (AGGR_MUL * x) : CoinMax (0.,(x+lb)/2));
      else                       return (direction ? AGGR_MUL : lb/AGGR_DIV);
    } else // [l,u]
      return (direction ? (x + (ub-x) / AGGR_DIV) : x - (x-lb) / AGGR_DIV);
  }
}


// Single fake tightening. Return
//
// -1   if infeasible
//  0   if no improvement
// +1   if improved
int fake_tighten (const CouenneCutGenerator *cg,
		  const OsiSolverInterface  *psi,
		  OsiCuts &cs, 
		  Bonmin::BabInfo * babInfo,

		  char direction,  // 0: left, 1: right
		  int index,       // index of the variable tested
		  const double *X, // point round which tightening is done
		  CouNumber *olb,  // cur. lower bound
		  CouNumber *oub,  //      upper
		  t_chg_bounds *chg_bds,
		  t_chg_bounds *f_chg) {
  int 
    ncols    = cg -> Problem () -> nVars (),
    objsense = cg -> Problem () -> Obj (0) -> Sense (),
    objind   = cg -> Problem () -> Obj (0) -> Body  () -> Index ();

  assert (objind >= 0);

  CouNumber 
    *lb       = cg -> Problem () -> Lb (),
    *ub       = cg -> Problem () -> Ub (),
    cutoff    = cg -> Problem () -> getCutOff (),
    xcur      = X [index],
    inner     = xcur,                                                 // point closest to current x
    //innerZ    = ((objsense == MINIMIZE) ? lb : ub) [objind],          // and associated dual bound
    outer     = (direction ? oub : olb) [index],                      // point closest to bound
    //outerZ    = COUENNE_INFINITY,                                     // ditto
    fb        = fictBounds (direction, xcur, lb [index], ub [index]); // starting point

  bool tightened = false;

  // This is a one-dimensional optimization problem between inner and
  // outer, on a monotone function of which we can compute the value
  // (with relative expense) but not the derivative.

#ifdef DEBUG
  CouNumber curdb     = ((objsense == MINIMIZE) ? lb : ub) [objind];  // current dual bound
  printf ("  x_%d.  x = %10g, lb = %g, cutoff = %g-----------------\n", index, xcur, curdb, cutoff);
#endif

  for (int iter = 0; iter < MAX_ITER; iter++) {

    (direction ? lb : ub) [index] = fb; 

#ifdef DEBUG
    char c1 = direction ? '-' : '>', c2 = direction ? '<' : '-';
    printf ("    #%3d: [%+10g -%c %+10g %c- %+10g] /\\/\\ ",iter,olb[index],c1,fb,c2, oub [index]);
    printf (" [%10g,%10g] <%g,%g>==> ",lb[index],ub[index],CoinMin(inner,outer),CoinMax(inner,outer));
#endif

    bool
      feasible  = btCore (cg, psi, cs, f_chg, babInfo, false), // true if feasible with fake bound
      betterbds = (objsense == MINIMIZE) ?                     // true if over cutoff
        (lb [objind] > cutoff) : 
        (ub [objind] < cutoff);

#ifdef DEBUG
    printf(" [%10g,%10g] lb = %g {fea=%d,btr=%d} ",lb[index],ub[index],lb[objind],feasible,betterbds);
#endif

    if (feasible && !betterbds) {

      // case 1: too tight, move inner out
      inner = fb;

      // restore initial bound
      CoinCopyN (chg_bds, ncols, f_chg);
      CoinCopyN (olb, ncols, lb);
      CoinCopyN (oub, ncols, ub);

    } else {

      // case 2: tightening valid, apply and move outer in

#ifdef DEBUG
      printf (" --> %cbound [x_%d]: %g --> %g",direction?'U':'L',index,(direction?oub:olb)[index],fb);
#endif

      // apply bound
      if (direction) {oub[index] = ub[index] = fb; chg_bds [index].setUpper (t_chg_bounds::CHANGED);}
      else           {olb[index] = lb[index] = fb; chg_bds [index].setLower (t_chg_bounds::CHANGED);}

      outer = fb; // we have at least a tightened bound, save it 

      tightened = true;

      // restore initial bound
      CoinCopyN (chg_bds, ncols, f_chg);
      CoinCopyN (olb, ncols, lb);
      CoinCopyN (oub, ncols, ub);

      //#if BR_TEST_LOG < 0 // for fair testing
      // check tightened problem for feasibility
      if (!(btCore (cg, psi, cs, chg_bds, babInfo, true))) {
#ifdef DEBUG
	printf ("\n    pruned by aggressive BT\n");
#endif
	return -1;
      }
      //#endif
    }

    // TODO: compute golden section
    fb = (inner + outer) / 2;

    //    if () fb = (          inner + (phi-1) * outer) / phi;
    //    else  fb = ((phi-1) * inner +           outer) / phi;

    //	if (!feasible)       
    //    fb = fictBounds (direction, xcur, 
    //	     direction ? lb [index] : outer,
    //	     direction ? outer      : ub [index]);
    //    fb = fictBounds (direction, xcur, CoinMin (inner, outer), CoinMax (inner, outer));

#ifdef DEBUG
    printf ("\n");
#endif
  }

#ifdef DEBUG
  printf ("\n");
  if (tightened) printf ("  [x%2d] pruned %s [%g, %g] -- lb = %g cutoff = %g\n", 
			 index,direction?"right":"left",olb[index], oub [index], lb [objind], cutoff);
#endif

  return tightened ? 1 : 0;
}

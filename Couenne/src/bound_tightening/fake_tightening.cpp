/*
 * Name:    fake_tightening.cpp
 * Author:  Pietro Belotti
 * Purpose: fake single bounds in variables to exclude parts of the solution space 
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneProblem.hpp"
#include "BonBabInfos.hpp"

//#define DEBUG

#define MAX_ITER  2 // max # fake tightening (inner) iterations 
#define AGGR_MUL  2 // the larger,  the more conservative. Must be > 0
#define AGGR_DIV  2 // the smaller, the more conservative. Must be > 1

// golden ratio, used to find the ideal bound
const CouNumber phi = 0.5 * (1. + sqrt (5.));

//	if (Jnlst()->ProduceOutput(Ipopt::J_VECTOR, J_BOUNDTIGHTENING)) {
//	  Jnlst()->Printf(Ipopt::J_VECTOR, J_BOUNDTIGHTENING,


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
int CouenneProblem::
fake_tighten (char direction,  // 0: left, 1: right
	      int index,       // index of the variable tested
	      const double *X, // point round which tightening is done
	      CouNumber *olb,  // cur. lower bound
	      CouNumber *oub,  //      upper
	      t_chg_bounds *chg_bds,
	      t_chg_bounds *f_chg) const {
  int
    ncols    = nVars (),
    //objsense = Obj (0) -> Sense (),
    objind   = Obj (0) -> Body  () -> Index ();

  assert (objind >= 0);

  bool 
    tightened = false,
    intvar    = variables_ [index] -> isInteger ();

  CouNumber 
    xcur      = X [index],
    inner     = xcur,                                                 // point closest to current x
    outer     = (direction ? oub : olb) [index],                      // point closest to bound
    fb        = fictBounds (direction, xcur, Lb (index), Ub (index)); // starting point

  // This is a one-dimensional optimization problem between inner and
  // outer, on a monotone function of which we can compute the value
  // (with relative expense) but not the derivative.

#ifdef DEBUG
  CouNumber curdb     = Lb (objind);// : Ub (objind);  // current dual bound
  printf ("  x_%d.  x = %10g, lb = %g, cutoff = %g-----------------\n", index,xcur,curdb,getCutOff());
#endif

  for (int iter = 0; iter < MAX_ITER; iter++) {

    if (intvar) {

      if (!direction) {inner = floor (inner); outer = ceil  (outer);}
      else            {inner = ceil  (inner); outer = floor (outer);}

      if ( direction && (inner > outer) ||
	  !direction && (inner < outer)) {

	// apply bound
	if (direction) {oub[index] = Ub (index) = fb; chg_bds[index].setUpper(t_chg_bounds::CHANGED);}
	else           {olb[index] = Lb (index) = fb; chg_bds[index].setLower(t_chg_bounds::CHANGED);}

	tightened = true;

	// restore initial bound
	CoinCopyN (chg_bds, ncols, f_chg);
	CoinCopyN (olb, ncols, Lb ());
	CoinCopyN (oub, ncols, Ub ());

	break;
      }

      if (direction  && ((fb < inner) || (fb > outer)) ||
	  !direction && ((fb > inner) || (fb < outer)))
	fb = 0.5 * (inner + outer);
    }

    if (direction) {Lb (index) = fb; f_chg [index].setLower (t_chg_bounds::CHANGED);} 
    else           {Ub (index) = fb; f_chg [index].setUpper (t_chg_bounds::CHANGED);}

    //    (direction ? lb_ : ub_) [index] = fb; 

#ifdef DEBUG
    char c1 = direction ? '-' : '>', c2 = direction ? '<' : '-';
    printf ("    #%3d: [%+10g -%c %+10g %c- %+10g] /\\/\\ ",iter,olb[index],c1,fb,c2, oub [index]);
    printf (" [%10g,%10g]<%g,%g>=> ",Lb (index),Ub (index),CoinMin(inner,outer),CoinMax(inner,outer));
#endif

    bool
      feasible  = btCore (f_chg),             // true if feasible with fake bound
      betterbds = Lb (objind) > getCutOff (); // true if over cutoff
      

#ifdef DEBUG
    printf(" [%10g,%10g] lb = %g {fea=%d,btr=%d} ",
	   Lb (index), Ub (index), Lb (objind),feasible,betterbds);
#endif

    if (feasible && !betterbds) {

      // case 1: too tight, move inner out
      inner = fb;

      // restore initial bound
      CoinCopyN (chg_bds, ncols, f_chg);
      CoinCopyN (olb, ncols, Lb ());
      CoinCopyN (oub, ncols, Ub ());

    } else {

      // case 2: tightening valid, apply and move outer in

#ifdef DEBUG
      printf (" --> %cbound [x_%d]: %g --> %g",direction?'U':'L',index,(direction?oub:olb)[index],fb);
      if (optimum_ && 
	  ((!direction &&
	    (optimum_ [index] >= olb [index]) && 
	    (optimum_ [index] <= Lb (index) - COUENNE_EPS)) ||
	   (direction &&
	    (optimum_ [index] <= oub [index]) && 
	    (optimum_ [index] >= Ub (index) + COUENNE_EPS)))) {
	printf ("fake tightening cuts out optimum: x%d=%g in [%g,%g] but not in [%g,%g]\n",
		index, olb [index], oub [index], Lb (index), Ub (index));
      }
#endif

      /*bool do_not_tighten = false;

      // check if cut best known solution
      if (optimum_) {
	if (direction) {
	  if ((oub [index] > optimum_ [index]) && 
	      (fb          < optimum_ [index])) {
	    printf ("aggressive bt cuts optimum ub %d: %g < %g < %g\n", 
		    index, fb, optimum_ [index], oub [index]);
	    do_not_tighten = true;
	  }
	} else {
	  if ((olb [index] < optimum_ [index]) && 
	      (fb          > optimum_ [index])) {
	    printf ("aggressive bt cuts optimum lb %d: %g < %g < %g\n", 
		    index, olb [index], optimum_ [index], fb);
	    do_not_tighten = true;
	  }
	}
	}*/

      //if (!do_not_tighten) {

	// apply bound
      if (direction) {oub[index]=Ub (index) = fb; chg_bds [index].setUpper (t_chg_bounds::CHANGED);}
      else           {olb[index]=Lb (index) = fb; chg_bds [index].setLower (t_chg_bounds::CHANGED);}

      outer = fb; // we have at least a tightened bound, save it 

      tightened = true;
	//}

      // restore initial bound
      CoinCopyN (chg_bds, ncols, f_chg);
      CoinCopyN (olb, ncols, Lb ());
      CoinCopyN (oub, ncols, Ub ());

      //#if BR_TEST_LOG < 0 // for fair testing
      // check tightened problem for feasibility
      if (!(btCore (chg_bds))) {
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


  Jnlst()->Printf(Ipopt::J_VECTOR, J_BOUNDTIGHTENING, "\n");
  if (tightened) 
    Jnlst()->Printf(Ipopt::J_VECTOR, J_BOUNDTIGHTENING, 
		    "  [x%2d] pruned %s [%g, %g] -- lb = %g cutoff = %g\n", 
		    index,direction?"right":"left",
		    olb[index],oub[index], Lb (objind), getCutOff ());

  return tightened ? 1 : 0;
}

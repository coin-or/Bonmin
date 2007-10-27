/*
 * Name:    impliedBounds-exprQuad.cpp
 * Author:  Pietro Belotti
 * Purpose: inferring bounds on independent variables of an exprQuad
 *          given bounds on the auxiliary variable
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <set>

#include "CoinHelperFunctions.hpp"
#include "exprQuad.hpp"

//#define DEBUG

// Compute bounds of exprQuad taking into account the unbounded variables
void computeQuadFiniteBound (const exprQuad *,
			     CouNumber &, CouNumber &, 
			     CouNumber *, CouNumber *,
			     int &, int &);


/// implied bound processing for quadratic form upon change in lower-
/// and/or upper bound of w, whose index is wind

bool exprQuad::impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg) {

  // return false;

#ifdef DEBUG
  printf ("################ implied bounds: [%g,%g], ", l [wind], u [wind]);
  print (); printf ("\n");
#endif

  // Nevermind about the nonlinear part stored in arglist_...

  // get all variables 

  std::set <int> 
    lin (index_,   index_   + nlterms_),
    qui (qindexI_, qindexI_ + nqterms_),
    quj (qindexJ_, qindexJ_ + nqterms_),
    indUnion;

  // CAUTION: this relies on the first version of bound expressions
  // for quadratic form, i.e. the sum of bounds of independent terms

  // find all indices appearing in the expression
  std::set_union (qui .begin (), qui .end (), 
		  quj .begin (), quj .end (),
		  std::inserter (indUnion, indUnion.begin ()));

  std::set_union (indUnion .begin (), indUnion .end (), 
		  lin      .begin (), lin      .end (),
		  std::inserter (indUnion, indUnion.begin ()));

  CouNumber qMin, qMax;

  int indInfLo = -1, indInfUp = -1;

  // get bounds for nonlinear part of the sum
  expression *lb, *ub;
  exprSum::getBounds (lb, ub);

  qMin = (*lb) ();
  qMax = (*ub) ();

  delete lb;
  delete ub;

  if (qMin < -COUENNE_INFINITY) indInfLo = -2;
  if (qMax >  COUENNE_INFINITY) indInfUp = -2;

  if ((indInfLo == -2) && 
      (indInfUp == -2)) 
    return false;
  
#ifdef DEBUG
  printf ("1st phase... inf=(%d,%d) q=[%g,%g].\n", indInfLo, indInfUp, qMin, qMax);
  for (std::set <int>:: iterator iter = indUnion.begin ();
       iter != indUnion.end (); ++iter) 
    printf ("%4d [%+6g %+6g]\n", *iter, l [*iter], u [*iter]);
#endif

  // compute bound on expression using only finite variable bounds

  computeQuadFiniteBound (this, qMin, qMax, l, u, indInfLo, indInfUp);

  // similar to impliedBounds-exprGroup. indInf* are -1 if no
  // variables are unbounded, i>0 if variable x_i is unbounded and -2
  // if at least two are, thus making implied bounds useless.

  qMin += c0_;
  qMax += c0_;

  if (((indInfLo == -2) && (indInfUp == -2)) || // at least two variables are unbounded
      ((indInfLo == -1) && (indInfUp == -1) &&  // or, none is (and qMin, qMax are exact) 
       (qMin > l [wind]) && (qMax < u [wind]))) // but explicit bounds are loose
    return false;

#ifdef DEBUG
  printf ("2nd phase... inf=(%d,%d) q=[%g,%g].\n", indInfLo, indInfUp, qMin, qMax);
#endif

  //////////////////////////////////////////////////////////////////////////////////
  //
  // now fill in b_i (constant term in both linear coefficient vectors)

  // prepare data structure for scanning all variables
  int
    minindex = *(indUnion.begin  ()), // minimum index
    maxindex = *(indUnion.rbegin ()), // maximum index
    nvars = maxindex - minindex + 1;  // upper bound on # variables involved

  CouNumber 
    *linCoeMin = new CouNumber [nvars], // min coeff of var x_i
    *linCoeMax = new CouNumber [nvars], // max coeff of var x_i
    *qii       = new CouNumber [nvars], // quadratic coeff
    *bCutLb    = new CouNumber [nvars], // amount to be drawn from the overall lower bound
    *bCutUb    = new CouNumber [nvars]; // amount to be drawn from the overall upper bound

  // fill all with 0
  CoinFillN (linCoeMin, nvars, 0.);
  CoinFillN (linCoeMax, nvars, 0.);
  CoinFillN (qii,       nvars, 0.);
  CoinFillN (bCutLb,    nvars, 0.);
  CoinFillN (bCutUb,    nvars, 0.);

  // assume all coefficients are finite
  for (int i=0; i<nlterms_; i++) {

    int ind = index_ [i];

    CouNumber
      coe = coeff_ [i], 
      li  = l [ind], 
      ui  = u [ind];

    ind -= minindex;

    linCoeMin [ind] += coe;
    linCoeMax [ind] += coe;

    if (coe > 0) { // contribution of linear term to bound depends on its coefficient
      if (li > -COUENNE_INFINITY) bCutLb [ind] += coe * li;
      if (ui <  COUENNE_INFINITY) bCutUb [ind] += coe * ui;
    } else {
      if (ui <  COUENNE_INFINITY) bCutLb [ind] += coe * ui;
      if (li > -COUENNE_INFINITY) bCutUb [ind] += coe * li;
    }
  }

#ifdef DEBUG
  printf ("linear filling (%d,%d): -----------------------\n", minindex, maxindex);
  for (std::set <int>:: iterator iter = indUnion.begin ();
       iter != indUnion.end (); ++iter) 
    printf ("%4d [%+6g %+6g] [%+6g %+6g]\n", *iter, 
	    linCoeMin [*iter - minindex], linCoeMax [*iter - minindex],
	    bCutLb    [*iter - minindex], bCutUb    [*iter - minindex]);
#endif

  // fill in remaining linear coefficients and quadratic ones
  for (int i=0; i<nqterms_; i++) {

    int qi = qindexI_ [i],
        qj = qindexJ_ [i];

    CouNumber coe = qcoeff_ [i],
      li = l [qi], lj = l [qj], 
      ui = u [qi], uj = u [qj];

    if (qi == qj) { // quadratic term

      qi -= minindex;

      qii [qi] = coe; // quadratic term

      CouNumber
	maxbUb = CoinMax (fabs (li), fabs (ui)),
	maxbLb = (li >= 0) ? (li) : (ui <= 0) ? (ui) : 0;

      if (maxbUb > COUENNE_INFINITY) maxbUb = 0;

      maxbUb *= maxbUb * coe;
      maxbLb *= maxbLb * coe;

      if (coe > 0) {
	bCutUb [qi] += maxbUb;
	bCutLb [qi] += maxbLb;
      } else {
	bCutUb [qi] += maxbLb;
	bCutLb [qi] += maxbUb;
      }
    } else { // product term

      coe *= 2;

      CouNumber *b1, *b2;

      if (coe > 0) {b1 = l; b2 = u;} 
      else         {b1 = u; b2 = l;}

      b1 += minindex;
      b2 += minindex;

      qi -= minindex;
      qj -= minindex;

      linCoeMin [qi] += coe * b1 [qj];
      linCoeMin [qj] += coe * b1 [qi];

      linCoeMax [qi] += coe * b2 [qj];
      linCoeMax [qj] += coe * b2 [qi];

      CouNumber
	addLo = CoinMin (CoinMin (li*lj, ui*uj),
		         CoinMin (ui*lj, li*uj)),
	addUp = CoinMax (CoinMax (li*lj, ui*uj), 
  	  	         CoinMax (ui*lj, li*uj));

      if (addLo < -COUENNE_INFINITY) addLo = 0;
      if (addUp >  COUENNE_INFINITY) addUp = 0;

      addLo *= coe;
      addUp *= coe;

      if (coe > 0) {
	bCutLb [qi] += addLo; bCutUb [qi] += addUp;
	bCutLb [qj] += addLo; bCutUb [qj] += addUp;
      } else {
	bCutLb [qi] += addUp; bCutUb [qi] += addLo;
	bCutLb [qj] += addUp; bCutUb [qj] += addLo;
      }
    }
  }

#ifdef DEBUG
  printf ("quad filling: -----------------------\n");
  for (std::set <int>:: iterator iter = indUnion.begin ();
       iter != indUnion.end (); ++iter) 
    printf ("%4d [%+6g %+6g] [%+6g %+6g]\n", *iter, 
	    linCoeMin [*iter - minindex], linCoeMax [*iter - minindex],
	    bCutLb    [*iter - minindex], bCutUb    [*iter - minindex]);
#endif

  // Done filling vectors /////////////////////////////////////////////////////////////
  // Now improve each independent variable in the set indUnion

  bool updated = false;

  for (std::set <int>:: iterator iter = indUnion.begin ();
       iter != indUnion.end (); ++iter) {

    int ind  = *iter, 
        indn = ind - minindex;

    CouNumber 
      al = linCoeMin [indn],
      au = linCoeMax [indn],
      q  = qii       [indn];

#ifdef DEBUG
    CouNumber
      ol = l [ind], 
      ou = u [ind];
#endif

    if (fabs (q) < COUENNE_EPS) { // case 1: qii is zero, term is "linear"

      if ((al > 0) || (au < 0)) { // otherwise, min and max lin coe
				  // contain zero, not much to do...

	//          l         <= b_i x_i + c <= u
	// <===>    l - c_MAX <= b_i x_i     <= u - c_MIN
	//
	// c_MAX = qMax - bCutUb [indn]
	// c_MIN = qMin - bCutLb [indn]

	CouNumber 
	  l_b = l [wind] - qMax + bCutUb [indn],
	  u_b = u [wind] - qMin + bCutLb [indn];

	if (al > 0) { // care about min -- see outer if, it means 0 < al < au

	  if ((indInfUp == -1) || (indInfUp == ind))
	    updated = updateBound (-1, l + ind, (l_b) / ((l_b < 0) ? al : au)) || updated;
	  if ((indInfLo == -1) || (indInfLo == ind)) 
	    updated = updateBound (+1, u + ind, (u_b) / ((u_b < 0) ? au : al)) || updated;

#ifdef DEBUG
	  if (l [ind] > ol) printf ("0. l%d: %g --> %g\n", ind, ol, l [ind]);
	  if (u [ind] < ou) printf ("0. u%d: %g --> %g\n", ind, ou, u [ind]);
#endif
	} else { // only look at max, as al < au < 0

	  if ((indInfLo == -1) || (indInfLo == ind)) 
	    updated = updateBound (-1, l + ind, (u_b) / ((u_b < 0) ? al : au)) || updated;
	  if ((indInfUp == -1) || (indInfUp == ind)) 
	    updated = updateBound (+1, u + ind, (l_b) / ((l_b < 0) ? au : al)) || updated;

#ifdef DEBUG
	  if (l [ind] > ol) printf ("1. l%d: %g --> %g\n", ind, ol, l [ind]);
	  if (u [ind] < ou) printf ("1. u%d: %g --> %g\n", ind, ou, u [ind]);
#endif
	}
      }

    } else if (q > 0) {

      // skip if constant term is unbounded from below

      if ((indInfLo != -1) && 
	  (indInfLo != ind)) 
	continue;

      // case 2: qii is positive, the parabola is facing upwards and
      // we only need to look at w_U = u [wind]

      // there are two parabolas, with linear coefficient linCoeMin
      // and linCoeMax, respectively. If both cut the line $w=w_U$
      // then take the minimum of the cut points as new lower bound
      // and similarly for the new upper bound.

      // Fortunately there are just two values of linCoe which we have
      // to look at, as the parabola is of the form $q_{ii} x_i^2 +
      // \hat b_i x_i + c \le w_u$ and the solution contains $\sqrt
      // {b_i^2 - 4q_{ii}(c-u)}$, whose maximum is attained with b_i
      // maximum or minimum, that is, at linCoeMax and linCoeMin.

      CouNumber 
	deltaSecond = 4 * q * (qMin - bCutLb [indn] - u [wind]),
	deltaUp     = au*au - deltaSecond,
	deltaLo     = al*al - deltaSecond;

      // First case, both parabolas cut the line

      if ((deltaUp >= 0) && 
	  (deltaLo >= 0)) {

	updated = updateBound (-1, l + ind, (- au - sqrt (deltaUp)) / (2*q)) || updated;
	updated = updateBound (+1, u + ind, (- al + sqrt (deltaLo)) / (2*q)) || updated;

#ifdef DEBUG
	if (l [ind] > ol) printf ("2. l%d: %g --> %g\n", ind, ol, l [ind]);
	if (u [ind] < ou) printf ("2. u%d: %g --> %g\n", ind, ou, u [ind]);
#endif

      } else if (deltaUp >= 0) { // only  left-leaning parabola does

	updated = updateBound (-1, l + ind, (- au - sqrt (deltaUp)) / (2*q)) || updated;
	updated = updateBound (+1, u + ind, (- au + sqrt (deltaUp)) / (2*q)) || updated;
	//	addCoeffCut ();

#ifdef DEBUG
	if (l [ind] > ol) printf ("3. l%d: %g --> %g\n", ind, ol, l [ind]);
	if (u [ind] < ou) printf ("3. u%d: %g --> %g\n", ind, ou, u [ind]);
#endif

      } else if (deltaLo >= 0) { // only right-leaning parabola does

	updated = updateBound (-1, l + ind, (- al - sqrt (deltaLo)) / (2*q)) || updated;
	updated = updateBound (+1, u + ind, (- al + sqrt (deltaLo)) / (2*q)) || updated;
	//	addCoeffCut ();

#ifdef DEBUG
	if (l [ind] > ol) printf ("4. l%d: %g --> %g\n", ind, ol, l [ind]);
	if (u [ind] < ou) printf ("4. u%d: %g --> %g\n", ind, ou, u [ind]);
#endif

      } else { // none of them does, the problem is infeasible

	updated = updateBound (-1, l+ind, +1) || updated;
	updated = updateBound (+1, u+ind, -1) || updated;

#ifdef DEBUG
	printf ("5. infeasible!\n");
#endif
      }

      // if only one parabola cuts the line, take its lower and upper
      // bounds as the new bounds

      // TODO: in the 2nd and 3rd cases, add constraint $\hat b \ge
      // \sqrt {4q_{ii}(c-u)}$ or the opposite, to ensure that
      // parabola is at least tangent to line w=w_U

      // Also, lower bound can be of help if x_i has a bound between
      // the two intersections of the parabola with line w=w_L

    } else {

      // case 3: qii is negative, the parabola is facing downwards and
      // we only need to look at l [wind]

      // skip if constant term is unbounded from above
      if ((indInfUp != -1) && 
	  (indInfUp != ind)) 
	continue;

      CouNumber 
	deltaSecond = 4 * q * (qMax - bCutUb [indn] - l [wind]),
	deltaUp     = au*au - deltaSecond,
	deltaLo     = al*al - deltaSecond;

      // First case, both parabolas cut the line

      if ((deltaUp >= 0) && 
	  (deltaLo >= 0)) {

	updated = updateBound (-1, l + ind, (al - sqrt (deltaLo)) / (-2*q)) || updated;
	updated = updateBound (+1, u + ind, (au + sqrt (deltaUp)) / (-2*q)) || updated;

#ifdef DEBUG
	if (l [ind] > ol) printf ("6. l%d: %g --> %g\n", ind, ol, l [ind]);
	if (u [ind] < ou) printf ("6. u%d: %g --> %g\n", ind, ou, u [ind]);
#endif

      } else if (deltaUp >= 0) { // only  left-leaning parabola does

	updated = updateBound (-1, l + ind, (au - sqrt (deltaUp)) / (-2*q)) || updated;
	updated = updateBound (+1, u + ind, (au + sqrt (deltaUp)) / (-2*q)) || updated;
	//	addCoeffCut ();

#ifdef DEBUG
	if (l [ind] > ol) printf ("7. l%d: %g --> %g\n", ind, ol, l [ind]);
	if (u [ind] < ou) printf ("7. u%d: %g --> %g\n", ind, ou, u [ind]);
#endif

      } else if (deltaLo >= 0) { // only right-leaning parabola does

	updated = updateBound (-1, l + ind, (al - sqrt (deltaLo)) / (-2*q)) || updated;
	updated = updateBound (+1, u + ind, (al + sqrt (deltaLo)) / (-2*q)) || updated;
	//	addCoeffCut ();
#ifdef DEBUG
	if (l [ind] > ol) printf ("8. l%d: %g --> %g\n", ind, ol, l [ind]);
	if (u [ind] < ou) printf ("8. u%d: %g --> %g\n", ind, ou, u [ind]);
#endif

      } else { // none of them does, the problem is infeasible

	updated = updateBound (-1, l+ind, +1) || updated;
	updated = updateBound (+1, u+ind, -1) || updated;

#ifdef DEBUG
	printf ("9. infeasible\n");
#endif
      }
    }
  }

  delete [] linCoeMin;
  delete [] linCoeMax;
  delete [] qii;
  delete [] bCutLb;
  delete [] bCutUb;

  return updated;
}

/*
 * Name:    decomposeTerm.cpp
 * Author:  Pietro Belotti
 * Purpose: decompose sums and products
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneProblem.hpp"
#include "exprMul.hpp"
#include "exprPow.hpp"
#include "exprQuad.hpp"
#include "lqelems.hpp"

//#define DEBUG

/// given (expression *) element of sum, returns (coe,ind0,ind1)
/// depending on element:
///
/// 1) a * x_i ^ 2   ---> (a,i,?)   return COU_EXPRPOW
/// 2) a * x_i       ---> (a,i,?)   return COU_EXPRVAR
/// 3) a * x_i * x_j ---> (a,i,j)   return COU_EXPRMUL
/// 4) a             ---> (a,?,?)   return COU_EXPRCONST
///
/// x_i and/or x_j may come from standardizing other (linear or
/// quadratic operator) sub-expressions

void CouenneProblem::decomposeTerm (expression *term,
				    CouNumber initCoe,
				    CouNumber &c0,
				    LinMap  &lmap,
				    QuadMap &qmap) {

  switch (term -> code ()) {

    // easy cases ////////////////////////////////////////////////////////////////////////

  case COU_EXPRCONST: /// a constant
    c0 += initCoe * term -> Value ();
    break;

  case COU_EXPRVAR:   /// a variable
    lmap.insert (term -> Index (), initCoe);
    break;

  case COU_EXPROPP:   /// the opposite of a term
    decomposeTerm (term -> Argument (), -initCoe, c0, lmap, qmap);
    break;

  case COU_EXPRSUB:   /// a subtraction
    decomposeTerm (term -> ArgList () [0],  initCoe, c0, lmap, qmap);
    decomposeTerm (term -> ArgList () [1], -initCoe, c0, lmap, qmap);
    break;

  case COU_EXPRQUAD: { /// a quadratic form

    exprQuad *t = dynamic_cast <exprQuad *> (term);
    exprQuad::sparseQ &M = t -> getQ ();

    for (exprQuad::sparseQ::iterator row = M.begin (); 
	 row != M.end (); ++row) {

      int xind = row -> first -> Index ();

      for (exprQuad::sparseQcol::iterator col = row -> second.begin (); 
	   col != row -> second.end (); ++col) {
	qmap.insert (xind, col -> first -> Index (), initCoe * col -> second);
      }
    }
  } // NO break here, exprQuad generalizes exprGroup

  case COU_EXPRGROUP: { /// a linear term

    exprGroup *t = dynamic_cast <exprGroup *> (term);
    exprGroup::lincoeff &lcoe = t -> lcoeff ();

    //  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el)
    for (int n = lcoe.size (), i=0; n--; i++)
      lmap.insert (lcoe [i].first -> Index (), initCoe * lcoe [i].second);

    c0 += initCoe * t -> getc0 ();
  } // NO break here, exprGroup generalizes exprSum

  case COU_EXPRSUM: { /// a sum of (possibly) nonlinear elements

    expression **al = term -> ArgList ();
    for (int i = term -> nArgs (); i--;)
      decomposeTerm (*al++, initCoe, c0, lmap, qmap);

  } break;

  // not-so-easy cases /////////////////////////////////////////////////////////////////
  //
  // cannot add terms as it may fill up the triplet

  case COU_EXPRMUL: { /// a product of n factors /////////////////////////////////////////

    std::map <int, CouNumber> indices;
    CouNumber coe = initCoe;

    // return list of variables (some of which auxiliary)
    flattenMul (term, coe, indices);

#ifdef DEBUG
    printf ("from flattenmul: [%g] ", coe);
    for (std::map <int, CouNumber>::iterator itt = indices.begin ();
	 itt != indices.end(); ++itt)
      printf (" %d,%g",
	      itt -> first,
	      itt -> second);
    printf ("\n");
#endif

    // based on number of factors, decide what to return
    switch (indices.size ()) {

    case 0: // no variables in multiplication (hmmm...)
      c0 += coe;
      break;

    case 1: { // only one term (may be with !=1 exponent)

      std::map <int, CouNumber>::iterator one = indices.begin ();
      int       index = one -> first;
      CouNumber expon = one -> second;

      if      (fabs (expon - 1) < COUENNE_EPS) lmap.insert (index, coe);
      else if (fabs (expon - 2) < COUENNE_EPS) qmap.insert (index, index, coe);
      else {
	exprAux *aux = addAuxiliary 
	  (new exprPow (new exprClone (Var (index)),
			new exprConst (expon)));

	//linsert (lmap, aux -> Index (), initCoe); // which of these three is correct?
	//linsert (lmap, aux -> Index (), initCoe * coe);
	lmap.insert (aux -> Index (), coe);
      }
    } break;

    case 2: { // two terms

      int ind0, ind1;

      std::map <int, CouNumber>::iterator one = indices.begin (), 
	two = one;
      ++two; // now "two" points to the other variable

      // first variable
      if (fabs (one -> second - 1) > COUENNE_EPS) {
	exprAux *aux = addAuxiliary (new exprPow (new exprClone (Var (one -> first)),
						  new exprConst (one -> second)));
	ind0 = aux -> Index ();
      } else ind0 = one -> first;

      // second variable
      if (fabs (two -> second - 1) > COUENNE_EPS) {
	exprAux *aux = addAuxiliary (new exprPow (new exprClone (Var (two -> first)),
						  new exprConst (two -> second)));
	ind1 = aux -> Index ();
      } else ind1 = two -> first;

      qmap.insert (ind0, ind1, coe);
    } break;

    default: { 

      // create new auxiliary variable containing product of 3+ factors

      expression **al = new expression * [indices.size ()];
      std::map <int, CouNumber>::iterator one = indices.begin ();

      for (int i=0; one != indices.end (); ++one, i++) 
	if (fabs (one -> second - 1) > COUENNE_EPS) {
	  exprAux *aux = addAuxiliary (new exprPow (new exprClone (Var (one -> first)),
						    new exprConst (one -> second)));
	  al [i] = new exprClone (aux);
	} else al [i] = new exprClone (Var (one -> first));

      // TODO: when we have a convexification for \prod_{i \in I}...
      //      exprAux *aux = addAuxiliary (new exprMul (al, indices.size ()));

      exprMul *mul = new exprMul (al, indices.size ());
      exprAux *aux = mul -> standardize (this);
      lmap.insert (aux -> Index (), coe);

    } break;
    }

  } break; // end of case COU_EXPRMUL

  case COU_EXPRPOW: { // expression = f(x)^g(x) ////////////////////////////////////////////////

    expression **al = term -> ArgList (); 

    if (al [1] -> Type () != CONST) { 

      // non-constant exponent, standardize the whole term and add
      // linear component (single aux)

      expression *aux = term -> standardize (this);
      if (!aux) aux = term;
      lmap.insert (aux -> Index (), initCoe);

    } else { // this is of the form f(x)^k.  If k=2, return square. If
	     // k=1, return var. Otherwise, generate new auxiliary.

      expression *aux = (*al) -> standardize (this);

      if (!aux)
	aux = *al; // it was a simple variable, and was not standardized.

      CouNumber expon = al [1] -> Value ();
      int ind = aux -> Index ();

      if      (fabs (expon - 1) < COUENNE_EPS) lmap.insert (ind, initCoe);
      else if (fabs (expon - 2) < COUENNE_EPS) qmap.insert (ind, ind, initCoe);
      else {
	exprAux *aux2 = addAuxiliary 
	  (new exprPow (new exprClone (aux), new exprConst (expon))); // TODO: FIX!
	lmap.insert (aux2 -> Index (), initCoe);
      }
    }
  } break;

  default: { /// otherwise, simply standardize expression 

    expression *aux = term -> standardize (this);
    if (!aux) 
      aux = term;
    lmap.insert (aux -> Index (), initCoe);
  } break;
  }
}

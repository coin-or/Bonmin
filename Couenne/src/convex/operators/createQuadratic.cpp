/*
 * Name:    createQuadratic.cpp
 * Author:  Pietro Belotti
 * Purpose: decompose sums and products
 *
 * This file is licensed under the Common Public License (CPL)
 */


#include <CouenneTypes.h>
#include <CouenneProblem.hpp>
#include <exprMul.hpp>
#include <exprPow.hpp>
#include <exprQuad.hpp>


/// re-organizes multiplication and stores indices (and exponents) of its variables
void flattenMul (expression *mul, CouNumber &coe, 
		 std::map <int, CouNumber> &indices, 
		 CouenneProblem *p) {

  int nargs = mul -> nArgs ();
  expression **al = mul -> ArgList ();

  // for each factor (variable, function, or constant) of the product
  for (int i=0; i < nargs; i++) { 

    expression *arg = al [i];

    switch (arg -> code ()) {

    case COU_EXPRCONST: // change scalar multiplier

      coe *= arg -> Value ();
      break;

    case COU_EXPRMUL:  // apply recursively

      flattenMul (arg, coe, indices, p);
      break;

    case COU_EXPRVAR: { // insert index or increment 

      std::map <int, CouNumber>::iterator 
	where = indices.find (arg -> Index ());

      if (where == indices.end ()) 
	indices.insert (std::pair <int, CouNumber> (arg -> Index (), 1));
      else ++ (where -> second);
    } break;

    case COU_EXPROPP: // equivalent to multiplying by -1

      coe *= -1;
      flattenMul (arg -> Argument (), coe, indices, p);
      break;

    case COU_EXPRPOW: {

      expression *base     = arg -> ArgList () [0],
	         *exponent = arg -> ArgList () [1];

      if (exponent -> Type () == CONST) { // could be of the form k x^2

	double expnum = exponent -> Value ();

	expression *aux = base -> standardize (p);

	if (!aux)
	  aux = base;

	std::map <int, CouNumber>::iterator 
	  where = indices.find (aux -> Index ());

	if (where == indices.end ())
	  indices.insert (std::pair <int, CouNumber> (aux -> Index (), expnum));
	else (where -> second += expnum);

	break;
      }  // otherwise, revert to default
    }

    default: { // for all other expression, add associated new auxiliary

      exprAux *aux = arg -> standardize (p);

      int ind = (aux) ? aux -> Index () : arg -> Index ();

      std::map <int, CouNumber>::iterator where = indices.find (ind);

      if (where == indices.end ()) 
	indices.insert (std::pair <int, CouNumber> (ind, 1));
      else ++ (where -> second);
    } break;
    }
  }
}


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

void decomposeTerm (CouenneProblem *p, expression *term,
		    CouNumber initCoe,
		    CouNumber &c0,
		    std::map <int,                 CouNumber> &lmap,
		    std::map <std::pair <int,int>, CouNumber> &qmap) {

  switch (term -> code ()) {

    // easy cases ////////////////////////////////////////////////////////////////////////

  case COU_EXPRCONST: /// a constant
    c0 += initCoe * term -> Value ();
    break;

  case COU_EXPRVAR:   /// a variable
    linsert (lmap, term -> Index (), initCoe);
    break;

  case COU_EXPROPP:   /// the opposite of a term
    decomposeTerm (p, term -> Argument (), -initCoe, c0, lmap, qmap);
    break;

  case COU_EXPRSUB:   /// a subtraction
    decomposeTerm (p, term -> ArgList () [0],  initCoe, c0, lmap, qmap);
    decomposeTerm (p, term -> ArgList () [1], -initCoe, c0, lmap, qmap);
    break;

  case COU_EXPRQUAD: { /// a quadratic form

    exprQuad *t = dynamic_cast <exprQuad *> (term);

    int       *qi = t -> getQIndexI ();
    int       *qj = t -> getQIndexJ ();
    CouNumber *qc = t -> getQCoeffs ();

    for (int i = t -> getnQTerms (); i--;)
      qinsert (qmap, *qi++, *qj++, initCoe * *qc++);
  } // NO break here, exprQuad generalizes exprGroup

  case COU_EXPRGROUP: { /// a linear term

    exprGroup *t = dynamic_cast <exprGroup *> (term);

    int       *ind = t -> getIndices ();
    CouNumber *coe = t -> getCoeffs ();

    for (int i = t -> getnLTerms (); i--;)
      linsert (lmap, *ind++, initCoe * *coe++);

    c0 += initCoe * t -> getc0 ();
  } // NO break here, exprGroup generalizes exprSum

  case COU_EXPRSUM: { /// a sum of (possibly) nonlinear elements

    expression **al = term -> ArgList ();
    for (int i = term -> nArgs (); i--;)
      decomposeTerm (p, *al++, initCoe, c0, lmap, qmap);

  } break;

    // not-so-easy cases /////////////////////////////////////////////////////////////////
    // cannot add terms as it may fill up the triplet

  case COU_EXPRMUL: { /// a product of n factors /////////////////////////////////////////

    std::map <int, CouNumber> indices;
    CouNumber coe = initCoe;

    // return list of variables (some of which auxiliary)
    flattenMul (term, coe, indices, p);

    // based on number of factors, decide what to return
    switch (indices.size ()) {

    case 0: // no variables in multiplication (hmmm...)
      c0 += coe;
      break;

    case 1: { // only one term (may be with >1 exponent)

      std::map <int, CouNumber>::iterator one = indices.begin ();
      int       index = one -> first;
      CouNumber expon = one -> second;

      if      (fabs (expon - 1) < COUENNE_EPS) linsert (lmap, index, coe);
      else if (fabs (expon - 2) < COUENNE_EPS) qinsert (qmap, index, index, coe);
      else {
	exprAux *aux = p -> addAuxiliary 
	  (new exprPow (new exprClone (p -> Var (index)),
			new exprConst (expon)));

	linsert (lmap, aux -> Index (), initCoe);
      }
    } break;

    case 2: { // two terms

      int ind0, ind1;

      std::map <int, CouNumber>::iterator one = indices.begin (), 
	two = one;
      ++two; // now "two" points to the other variable

      // first variable
      if (fabs (one -> second - 1) > COUENNE_EPS) {
	exprAux *aux = p -> addAuxiliary (new exprPow (new exprClone (p -> Var (one -> first)),
						       new exprConst (one -> second)));
	ind0 = aux -> Index ();
      } else ind0 = one -> first;

      // second variable
      if (fabs (two -> second - 1) > COUENNE_EPS) {
	exprAux *aux = p -> addAuxiliary (new exprPow (new exprClone (p -> Var (two -> first)),
						       new exprConst (two -> second)));
	ind1 = aux -> Index ();
      } else ind1 = two -> first;

      qinsert (qmap, ind0, ind1, coe);
    } break;

    default: { 

      // create new auxiliary variable containing product of 3+ factors

      expression **al = new expression * [indices.size ()];
      std::map <int, CouNumber>::iterator one = indices.begin ();

      for (int i=0; one != indices.end (); one++, i++) 
	if (fabs (one -> second - 1) > COUENNE_EPS) {
	  exprAux *aux = p -> addAuxiliary (new exprPow (new exprClone (p -> Var (one -> first)),
							 new exprConst (one -> second)));
	  al [i] = new exprClone (aux);
	} else al [i] = new exprClone (p -> Var (one -> first));

      // TODO: when we have a convexification for \prod_{i \in I}...
      //      exprAux *aux = p -> addAuxiliary (new exprMul (al, indices.size ()));

      exprMul *mul = new exprMul (al, indices.size ());
      exprAux *aux = mul -> standardize (p);

      linsert (lmap, aux -> Index (), coe);

    } break;
    }

  } break; // end of case COU_EXPRMUL

  case COU_EXPRPOW: { // expression = f(x)^g(x) ////////////////////////////////////////////////


    expression **al = term -> ArgList (); 

    if (al [1] -> Type () != CONST) { 

      // non-constant exponent, standardize the whole term and add
      // linear component (single aux)

      expression *aux = term -> standardize (p);
      if (!aux) aux = term;
      linsert (lmap, aux -> Index (), initCoe);

    } else { // this is of the form f(x)^k.  If k=2, return square. If
	     // k=1, return var. Otherwise, generate new auxiliary.

      expression *aux = (*al) -> standardize (p);

      if (!aux)
	aux = *al; // it was a simple variable, and was not standardized.

      CouNumber expon = al [1] -> Value ();
      int ind = aux -> Index ();

      if      (fabs (expon - 1) < COUENNE_EPS) linsert (lmap, ind, initCoe);
      else if (fabs (expon - 2) < COUENNE_EPS) qinsert (qmap, ind, ind, initCoe);
      else {
	exprAux *aux = p -> addAuxiliary (term);
	linsert (lmap, aux -> Index (), initCoe);
      }
    }
  } break;

  default: { /// otherwise, simply standardize expression 

    expression *aux = term -> standardize (p);
    if (!aux) 
      aux = term;
    linsert (lmap, aux -> Index (), initCoe);
  } break;
  }
}

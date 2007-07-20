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

      std::map <int, CouNumber>::iterator where = indices.find (arg -> Index ());

      if (where == indices.end ()) 
	indices.insert (std::pair <int, CouNumber> (arg -> Index (), 1));
      else ++ (where -> second);
    } break;

    default: { // for all other expression, add associated new auxiliary

      exprAux *aux = p -> addAuxiliary (arg);

      std::map <int, CouNumber>::iterator where = indices.find (aux -> Index ());

      if (where == indices.end ()) 
	indices.insert (std::pair <int, CouNumber> (aux -> Index (), 1));
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

int decomposeTerm (CouenneProblem *p, expression *term, CouNumber &coe, int &ind0, int &ind1) {

  switch (term -> code ()) {

    // easy cases ////////////////////////////////////////////////////////////////////////
    // only add term to triplet

  case COU_EXPRCONST: /// a constant
    coe *= term -> Value ();
    return COU_EXPRCONST;

  case COU_EXPRVAR:   /// a variable
    if (ind0 < 0) {ind0 = term -> Index (); return COU_EXPRVAR;}
    else          {ind1 = term -> Index (); return (ind1 == ind0) ? COU_EXPRPOW : COU_EXPRMUL;}

  case COU_EXPROPP: ///
    coe = -coe;
    return decomposeTerm (p, term -> Argument (), coe, ind0, ind1);

    // not-so-easy cases /////////////////////////////////////////////////////////////////
    // cannot add terms as it may fill up the triplet

  case COU_EXPRMUL: { /// will probably have to use recursion

    std::map <int, CouNumber> indices;

    // return list of variables (some of which auxiliary)
    flattenMul (term, coe, indices, p);

    // based on number of factors, decide what to return
    switch (indices.size ()) {

    case 0: return COU_EXPRCONST; // no variables in multiplication (hmmm...)

    case 1: { // only one term (may be with >1 exponent)

      std::map <int, CouNumber>::iterator one = indices.begin ();

      if (fabs (one -> second - 1) > COUENNE_EPS) ind0 = one -> first;
      else {
	exprAux *aux = p -> addAuxiliary (new exprPow (new exprClone (p -> Var (one -> first)),
						       new exprConst (one -> second)));
	ind0 = aux -> Index ();
      }

      return COU_EXPRVAR;
    }

    case 2: { // two terms

      std::map <int, CouNumber>::iterator one = indices.begin (), two = one;
      ++two;

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

      return (ind0 == ind1) ? COU_EXPRPOW : COU_EXPRMUL;
    }

    default: {// create new auxiliary variable containing product of 3+ factors

      expression **al = new expression * [indices.size ()];
      std::map <int, CouNumber>::iterator one = indices.begin ();

      for (int i=0; one != indices.end (); one++, i++) 
	if (fabs (one -> second - 1) > COUENNE_EPS) {
	  exprAux *aux = p -> addAuxiliary (new exprPow (new exprClone (p -> Var (one -> first)),
							 new exprConst (one -> second)));
	  al [i] = new exprClone (aux);
	} else al [i] = new exprClone (p -> Var (one -> first));

      exprAux *aux = p -> addAuxiliary (new exprMul (al, indices.size ()));

      ind0 = aux -> Index ();
      return COU_EXPRVAR;
    }
    }
  }      

  case COU_EXPRPOW: { // this is something of the form x^k.  If k=2,
		      // return square. If k=1, return var. Otherwise,
		      // generate new auxiliary.

    expression **al  = term -> ArgList (); 
    expression  *aux = (*al) -> standardize (p);

    if (!aux)
      aux = *al; // it was a simple variable, and was not standardized.

    // special case: exponent = 1
    if ((al [1] -> Type () == CONST) && 
	(fabs (al [1] -> Value () - 1) < COUENNE_EPS)) { // trivial power, 

      if (ind0 < 0) {ind0 = aux -> Index (); return COU_EXPRVAR;}
      else          {ind1 = aux -> Index (); return (ind1 == ind0) ? COU_EXPRPOW : COU_EXPRMUL;}
    }

    // behave differently if one space is occupied already

    if (ind0 >= 0) {

      // standardize power and put new aux as second argument
      //expression *aux = p -> addAuxiliary (term);
      //if (!aux) aux = term;

      ind1 = aux -> Index ();
      return (ind1 == ind0) ? COU_EXPRPOW : COU_EXPRMUL;

    } else {

      // special subcase, exponent is 2
      if ((al [1] -> Type () == CONST) && 
	  (fabs (al [1] -> Value () - 2) < COUENNE_EPS)) {

	//exprAux *aux = p -> addAuxiliary (*al);
	//if (!aux) 
	//aux = *al;

	ind0 = aux -> Index ();
	return COU_EXPRPOW;
      } else {

	//exprAux *aux = p -> addAuxiliary (term);
	//if (!aux)
	//  aux = term;

	ind0 = aux -> Index ();

	// standardize argument and return power
	return COU_EXPRVAR;
      }
    }
  }

  default: { /// for all other cases, standardize this expression and
	     /// return its data
    expression *aux = p -> addAuxiliary (term);
    if (!aux) 
      aux = term;

    if (ind0 >= 0) {ind1 = aux -> Index (); return (ind1 == ind0) ? COU_EXPRPOW : COU_EXPRMUL;}
    else           {ind0 = aux -> Index (); return COU_EXPRVAR;}    
  }
  }
}

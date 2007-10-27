/*
 * Name:    flattenMul.cpp
 * Author:  Pietro Belotti
 * Purpose: flatten multiplication expression tree into monomial
 *          c*\Prod_{k\in K} x_{i_k}^{p_k}
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneProblem.hpp"

/// re-organizes multiplication and stores indices (and exponents) of
/// its variables
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

      coe = -coe;
      flattenMul (arg -> Argument (), coe, indices, p);
      break;

    case COU_EXPRPOW: { // power

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

    case COU_EXPRSUM: // well, only if there is one element

      if (arg -> nArgs () == 1) {
	flattenMul (arg, coe, indices, p);
	break;
      }

    default: { // for all other expression, add associated new auxiliary

      exprAux *aux = arg -> standardize (p);

      int ind = (aux) ? aux -> Index () : arg -> Index ();

      std::map <int, CouNumber>::iterator 
	where = indices.find (ind);

      if (where == indices.end ()) 
	indices.insert (std::pair <int, CouNumber> (ind, 1));
      else ++ (where -> second);

    } break;
    }
  }
}

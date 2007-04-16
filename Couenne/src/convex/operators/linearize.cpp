/*
 * Name:    linearize.cpp
 * Author:  Pietro Belotti
 * Purpose: standardize sum expressions by flattening arguments
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprSum.h>
#include <exprGroup.h>
#include <exprConst.h>

#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>


/// used to create linear term
typedef struct {
  int index;
  CouNumber coeff;
} monomial;


/// flatten single expression, pass its coefficient/index to vectors
void flatten (expression *&arg, CouenneProblem *p, 
	      std::vector <monomial> *terms, 
	      CouNumber &a0, CouNumber sign) {

  // replace subexpression with its auxiliary
  register exprVar *subst;
  if ((subst = arg -> standardize (p)))
    arg = new exprClone (subst);

  // add corresponding pair (index,coeff) to the vector
  register int ind = arg -> Index ();

  if (ind >= 0) { // term is a non-constant monomial
    monomial term = {ind, sign};
    terms -> push_back (term);
  } 
  else a0 += sign * arg -> Value (); // term is constant
}


/// get constant multiplicator and aux index from a multiplication
/*void extractIndCoe (expression *, int &, CouNumber &) {

}*/


/// used to sort term vector
inline bool termSortPred (const monomial& lhs, const monomial& rhs)
{return lhs.index < rhs.index;}


// Create standard formulation of sums of expressions

exprAux *exprSum::standardize (CouenneProblem *p) {

  // first of all, standardize all operands
  //  exprOp::standardize (p);

  // now simply return NULL, (the caller will assume there is nothing
  // to change), as a sum is already standard
  //  return p -> addAuxiliary (this);

  // rather than flattening, another smarter way (and it prevents
  // DuplicateIndex exceptions) is to return this as a linear term.
  // This involves checking all operands of this sum to see if they
  // are of the type k*f(x), meaning each will turn into something
  // like k*w. In the end, all terms k1*w and k2*w should be grouped
  // as (k1+k2)*w.

  std::vector <monomial> terms;
  CouNumber a0 = 0;

  //printf ("Here we are:");
  //print (std::cout);
  //printf ("\n");

  // one by one standardize each argument of arglist_ into a pair
  // (index, coeff)

  for (int i = 0; i < nargs_; i++) {

    /// current argument
    expression *arg = arglist_ [i];

    //printf ("now curing: ");
    //arg -> print (std::cout);
    //printf ("\n");

    switch (arg -> code ()) {

    case COU_EXPRVAR:   // this is a variable, just pick the index
      {
	monomial term = {arg -> Index (), 1.};
	terms. push_back (term);
      }
      break;

    case COU_EXPRCONST: // a constant, add it to scalar term
      a0 += arg -> Value ();
      break;

      // TODO:

      /*case COU_EXPRMUL:   // multiplication k*f1(x)*h*f2(x)*f3(x)*...
      extractIndCoe (arglist_ [i], term.index, term.coeff); // separate k*h from f1(x)*f2(x)*f3(x)*...
      terms. push_back (term);
      break;*/

    case COU_EXPRGROUP: // just enqueue all linear term and wait for COU_EXPRSUM
      {
	exprGroup *e    = dynamic_cast <exprGroup *> (arg);
	int       *ind  = e -> getIndices ();
	CouNumber *coe  = e -> getCoeffs  ();

	while (*ind>=0) {
	  monomial term = {*ind++, *coe++};
	  terms. push_back (term);
	}

	a0 += e -> getc0 ();
      }

    case COU_EXPRSUM:   // decompose each argument of inner sum
      {
	expression **args = arglist_ [i] -> ArgList ();
	for (int j = arglist_ [i] -> nArgs (); j--;)
	  flatten (args [j], p, &terms, a0, +1.);
      }
      break;

    case COU_EXPRSUB:   // decompose both arguments of inner difference
      flatten (arglist_ [i] -> ArgList () [0], p, &terms, a0, +1.);
      flatten (arglist_ [i] -> ArgList () [1], p, &terms, a0, -1.);
      break;

    case COU_EXPROPP:   // decompose inner opposite
      flatten (*(arglist_ [i] -> ArgPtr ()), p, &terms, a0, -1.);
      break;

      // TODO:

      //case COU_EXPRDIV:   // three cases: x/y, k/y, x/k
      //  if (arglist_ [0] -> Type () <= CONST)

    default: 
      flatten (arglist_ [i], p, &terms, a0, +1.);
      break;
    }

    // useless now
    delete arglist_ [i]; arglist_ [i] = NULL;
  }

  // Take care of the exprGroup, if we are.
  if (code () == COU_EXPRGROUP) {

    exprGroup *e = dynamic_cast <exprGroup *> (this);

    int       *ind  = e -> getIndices ();
    CouNumber *coe  = e -> getCoeffs  ();
    a0             += e -> getc0 ();

    while (*ind >= 0) {
      monomial term = {*ind++, *coe++};
      terms.push_back (term);
    }
  }

  // Group all terms with the same variable, and create an exprGroup

  int nterms = terms.size ();

  int       *indices = (int       *) malloc (nterms * sizeof (int));  
  CouNumber *coeffs  = (CouNumber *) malloc (nterms * sizeof (CouNumber));

  // sort terms by index, to eliminate duplicates
  std::sort (terms.begin(), terms.end(), termSortPred);

  int prevind = -2, ndiff = 0;

  for (int i=0; i<nterms; i++) {

    //printf ("term %d/%d: %d, %e\n", i, nterms, 
    //    terms [i].index, terms [i].coeff);

    int ind = terms [i].index;

    if (ind == prevind)
      coeffs [ndiff] += terms [i].coeff;
    else {
      indices [ndiff]   = prevind = ind;
      coeffs  [ndiff++] = terms [i].coeff;
    }
  }

  indices = (int       *) realloc (indices, (1+ndiff) * sizeof (int));
  coeffs  = (CouNumber *) realloc (coeffs,     ndiff  * sizeof (CouNumber));

  indices [ndiff] = -1;

  expression **zero = new expression * [1];
  *zero = new exprConst (0.);

  return (p -> addAuxiliary (new exprGroup (a0, indices, coeffs, zero, 1)));
}

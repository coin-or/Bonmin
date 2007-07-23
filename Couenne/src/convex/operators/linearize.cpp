/*
 * Name:    linearize.cpp
 * Author:  Pietro Belotti
 * Purpose: standardize sum expressions by flattening arguments
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <exprSum.hpp>
#include <exprMul.hpp>
#include <exprGroup.hpp>
#include <exprQuad.hpp>
#include <exprConst.hpp>

#include <CouenneProblem.hpp>
#include <CouenneCutGenerator.hpp>

#if 0

/// used to create linear term
typedef struct {
  int index;
  CouNumber coeff;
} monomial;


/// get constant multiplicator and aux index from a multiplication
void extractIndCoe (CouenneProblem *p, expression *&mul, int &index, CouNumber &coeff) {

  // standardize expression
  exprAux *w;

  if ((w = mul -> standardize (p))) {
    //    delete mul;
    mul = new exprClone (w);
  } else { // it is a constant or a variable itself, but it cannot happen

    printf ("Couenne: strange, multiplication is actually a variable\n");
    index = mul -> Index ();
    coeff = 1;
    return;
  }

  // Assume mul's image is a product of two elements: if both are 
  // nonconstant, just return a pair with mul's index

  /*printf ("standardized mul into: ");
  w -> print (std::cout);  printf (" := ");
  w -> Image () -> print (std::cout);  printf ("\n");*/

  expression **alist = w -> Image () -> ArgList ();

  if (alist) {

    // first case: w = x*y. Just store w's index

    if ((alist [0] -> Type () != CONST) &&
	(alist [1] -> Type () != CONST)) {
      index = w -> Index ();
      coeff = 1;
      return;
    }

    // remind auxset that this is actually not used, at least here
    w -> decreaseMult (); 

    // look for coefficient. Standardized expression has only two
    // arguments, so it can be arglist_ [0] or arglist_ [1] -- or both

    // three cases: c1*y, x*c2, or c1*c2

    coeff = 1;
    index = alist [0] -> Index ();

    if (index < 0) { // first element is a constant
      index  = alist [1] -> Index (); 
      coeff *= alist [0] -> Value ();
    } else 
      coeff *= alist [1] -> Value ();

    if (index < 0) // second element is a constant
      coeff *= alist [1] -> Value ();
  }
  else if (mul -> Argument ()) {
    index = mul -> Argument () -> Index ();
    coeff = (index < 0) ? mul -> Argument () -> Value () : 1;
  } else {
    index = -1;
    coeff = 0;
  }
}


/// flatten single expression, pass its coefficient/index to vectors
void flatten (expression *&arg, CouenneProblem *p, 
	      std::vector <monomial> *terms, 
	      CouNumber &a0, CouNumber sign) {

  if (arg -> code () == COU_EXPRMUL) {

    monomial term = {-1, 0.};
    extractIndCoe (p, arg, term.index, term.coeff); 

    // now we have a variable (constant) and a coefficient (its value).
    if (sign<0) term.coeff = - term.coeff; // invert if it appears with a minus sign

    if (term.index < 0) a0 += term.coeff;  // add to constant if it's a constant
    else terms -> push_back (term);        // add monomial otherwise

    return;
  }

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


/// used to sort term vector
inline bool termSortPred (const monomial& lhs, const monomial& rhs)
{return lhs.index < rhs.index;}


// Create standard formulation of sums of expressions

exprAux *exprSum::standardize (CouenneProblem *p) {

  // rather than flattening, a better way (preventing DuplicateIndex
  // exceptions) is to return this as a linear term.  This involves
  // checking all operands of this sum to see if they are of the type
  // k*f(x), meaning each will turn into something like k*w. In the
  // end, all terms k1*w and k2*w should be grouped as (k1+k2)*w.

  std::vector <monomial> terms;
  CouNumber a0 = 0;

  if (nargs_ == 1) { // check special cases

    if (code () == COU_EXPRSUM) {

      // only one argument, just return its standardized form as
      // actual expression
      exprAux *ret = (*arglist_) -> standardize (p);
      if (ret) {
	*arglist_ = NULL;
	return ret;
      }
    } else if (code () == COU_EXPRGROUP) {

      exprGroup *eg = dynamic_cast <exprGroup *> (this);

      // if only one linear argument and no extra (that is, nl = 0 and
      // c0 = 0) then no standardization
      if ((fabs (eg -> getc0 ()) < COUENNE_EPS) && 
	  ((*arglist_) -> Linearity () <= ZERO)) {

	int *index = eg -> getIndices ();
	if ((*index >= 0) && (index [1] == -1)) { // bingo! It's of the form (0) + cx

	  expression *var = NULL;

	  for (std::vector <exprVar *>::iterator i = p -> Variables ().begin (); 
	       i != p -> Variables ().end (); i++)
	    if ((*i) -> Index () == *index)
	      var = new exprClone (*i);

	  /*for (std::vector <exprAux *>::iterator i = p -> Auxiliaries ().begin (); 
	       i != p -> Auxiliaries ().end (); i++)
	    if ((*i) -> Index () == *index)
	    var = new exprClone (*i);*/

	  // is variable indexed by *index original or auxiliary?
	  CouNumber coe = *(eg -> getCoeffs ());
	  if (fabs (coe-1) < COUENNE_EPS)
	    return    NULL;//p -> addAuxiliary (var);
	  else return p -> addAuxiliary (new exprMul (new exprConst (coe), var));
	}
      }
    }

    // don't check for quadratic terms, as there would be at most one
  }

  // check if this is a potential quadratic form, and return it if so
  exprAux *q = createQuadratic (p);
  if (q) return q;

  // one by one standardize each argument of arglist_ into a pair
  // (index, coeff)

  for (int i = 0; i < nargs_; i++) {

    /// current argument
    expression *arg = arglist_ [i];

    /*printf ("now curing: "); fflush (stdout);
    arg -> print (std::cout);
    printf ("\n");*/

    switch (arg -> code ()) {

    case COU_EXPRVAR:   // this is a variable, just pick the index
      {
	monomial term = {arg -> Index (), 1.};
	terms. push_back (term);
      }
      break;

    case COU_EXPRCONST: // a constant, add it to scalar term
      //printf ("found constant\n");
      a0 += arg -> Value ();
      break;

    case COU_EXPRMUL: // multiplication k*f1(x)*h*f2(x)*f3(x)*...
      {
	monomial term = {-1, 0.};
	extractIndCoe (p, arglist_ [i], term.index, term.coeff); 
	if (term.index < 0) a0 += term.coeff;
	else terms. push_back (term);
      }
      break;

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
      // no break, as an exprGroup is a derived class of what's below

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
      //printf ("it's something else\n");
      flatten (arglist_ [i], p, &terms, a0, +1.);
      break;
    }

    // useless now
    // delete arglist_ [i]; 
    //arglist_ [i] = NULL;
  }

  // Take care of the exprGroup
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

  /////////////////////////////////////////////////////////////////////////////////

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

  /*for (int i=0; i<nterms; i++) 
    printf ("<%d,%g> ", indices [i], coeffs [i]);

  printf ("[%d terms]\n", nterms); */

  indices = (int       *) realloc (indices, (1+ndiff) * sizeof (int));
  coeffs  = (CouNumber *) realloc (coeffs,     ndiff  * sizeof (CouNumber));

  indices [ndiff] = -1;

  expression **zero = new expression * [1];
  *zero = new exprConst (0.);

  exprGroup *eg = new exprGroup (a0, indices, coeffs, zero, 1);

  /*printf ("eg -> "); fflush (stdout);
    eg -> print (); printf ("\n");*/

  exprAux *ret = p -> addAuxiliary (eg);

  /*if (!ret) printf ("returning null\n");
  else {
    printf ("returning "); fflush (stdout);
    ret -> print (); printf (" := "); fflush (stdout); 
    ret -> Image () -> print (); printf ("\n");
    }*/

  return ret;
}

#endif

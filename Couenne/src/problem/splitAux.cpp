/*
 * Name:    splitAux.cpp
 * Author:  Pietro Belotti
 * Purpose: extract auxiliary variable from implicit equality constraint
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneProblemElem.hpp>
#include <CouenneProblem.hpp>

#include <exprSum.hpp>
#include <exprMul.hpp>
#include <exprGroup.hpp>

//#define DEBUG

/// given an element of a sum, check if it is a variable (possibly
/// with a coefficient) and return its index (and the coefficient) if
/// it has not been spotted as an auxiliary

void elementBreak (expression *, int &, CouNumber &);


/// split a constraint w - f(x) = c into w's index (it is returned)
/// and rest = f(x) + c

int splitAux (CouenneProblem *p, CouNumber rhs, 
	      expression *body, expression *&rest, bool *wentAux) {

  int auxInd = -1,          // index of the auxiliary to be extracted
    code = body -> code (); // type of expression

  expression **alist = body -> ArgList ();

  switch (code) { // constraint h(x) = 0 may be in different forms:
		  // subtraction, sum, linear group

    /////////////////////////////////////////////////////////////////////////////

  case COU_EXPRSUB: { 

    // simplest case, we have w-f(x)=rhs or f(x)-w=rhs or f(x)-g(x)=rhs

    int pos = 0;
    CouNumber coeff;

    auxInd = (*alist) -> Index ();

    elementBreak (*alist, auxInd, coeff); // check first element 

    if ((auxInd < 0) || 
	wentAux [auxInd] || 
	(body -> dependsOn (&auxInd, 1) > 1)) {

      elementBreak (alist [1], auxInd, coeff); // check second element

      if ((auxInd < 0) ||                       // no index found
	  wentAux [auxInd] ||                   // or, found but invalid
	  (body -> dependsOn (&auxInd, 1) > 1)) // or, variable depends upon itself
	return -1;

      //      printf ("no dependence0 %d\n", auxInd);

      pos = 1;
    }

    //////////////////////

    expression *clone = alist [1 - pos] -> clone (); // what remains is the "independent" expression

    expression *auxdef = 
      (fabs (coeff - 1) < COUENNE_EPS) ?          // if coefficient is 1
      (clone) :                                   //     do nothing
      new exprMul (new exprConst (coeff), clone); //     else multiply it by coefficient

    rest = (fabs (rhs) < COUENNE_EPS) ?                              // no extra constant?
      (auxdef) :                                                     // just put other argument
      (new exprSum (auxdef, new exprConst ((pos==1) ? -rhs : rhs))); // otherwise sum it with \pm rhs

  } break;

  ////////////////////////////////////////////////////////////////////////////

  //  case COU_EXPRQUAD: // TODO
  case COU_EXPRGROUP:
  case COU_EXPRSUM: {

    // an exprGroup or exprSum. Decompose the expression and
    // re-assemble (to look for right variable)

    // data structure to be used below if there is a linear term.
    // which specifies position within arrays (negative for linear
    // part of exprGroup, positive for all elements of exprSum)
    int maxindex = -1, *linind = NULL, nlin = 0, which = 1;
    CouNumber c0 = 0., *lincoe = NULL, auxcoe = 1;
    bool which_was_set = false;

    if (code != COU_EXPRSUM) { // check indices of linear part

      exprGroup *egBody = dynamic_cast <exprGroup *> (body);

      // import exprGroup linear structure

      linind = egBody -> getIndices ();
      lincoe = egBody -> getCoeffs  ();
      c0     = egBody -> getc0      ();

      for (int j; (j = linind [nlin]) >= 0; nlin++)
	if ((j > maxindex) && 
	    !(wentAux [j]) && 
	    (fabs (lincoe [nlin]) > COUENNE_EPS) &&
	    (body -> dependsOn (p, &j, 1) <= 1)) {

	  //	  printf ("no dependence (1) %d on", j);
	  //	  body -> print (); printf ("\n");

	  which    = - nlin - 1;    // mark which with negative number 
	  auxcoe   = lincoe [nlin];
	  maxindex = j;
	}
    }

    if (which != 1) 
      which_was_set = true;
    else which = -1;

    // check indices of possible linear elements of (nonlinear) sum

    for (int i = body -> nArgs (); i--;) {

      CouNumber coeff;
      int index;

      elementBreak (alist [i], index, coeff);

      if ((index > maxindex) &&
	  !(wentAux [index]) &&
	  (fabs (coeff) > COUENNE_EPS) && 
	  (body -> dependsOn (&index, 1) <= 1)) {

	//	printf ("no dependence2 %d\n", index);

	maxindex = index;
	which    = i;
	auxcoe   = coeff;
      }
    }

    if (!which_was_set && (which == -1)) // which has not been set
      return -1;

    ///////////////////////////////////////////////////////////////////

    if (maxindex < 0) break; // no substitute found, bail out of this expression

    // create a new exprGroup or exprSum with all elements but the
    // extracted auxiliary

    // start with exprSum
    int nargs = body -> nArgs ();
    expression **newarglist;

#ifdef DEBUG
    printf (" [[ind %d, coe %.1g, wh %d, nargs %d, nlin %d]] ", 
	    maxindex, auxcoe, which, nargs, nlin);
#endif

    if (nargs > 0) { // there is an element in the nonlinear sum to be drawn

      int j, mid = (which < 0) ? nargs : which;

      newarglist = new expression * [nargs + 1];

      for (j=0; j<mid;   j++) newarglist [j]   = alist [j];// newarglist [j]->  print ();printf(",");}
      for (j++; j<nargs; j++) newarglist [j-1] = alist [j];// newarglist [j-1]->print ();printf(";");}

      // nl arglist is done, later decide whether to incorporate it as
      // it is or with a coefficient

    } else { // no nonlinear arguments, or the only one is the new aux
      nargs++; // !!!!!!!!!!!!!!!!!!!!!!!!!
      newarglist  = new expression *;
      *newarglist = new exprConst (0);
    }

    // form rhs linear part ////////////////////////////////////////////////////

    int       *linind2 = NULL;
    CouNumber *lincoe2 = NULL;

    if (nlin > 0) { // there is an element in the linear sum to be drawn

      int mid = (which >= 0) ? nlin : - which - 1;

      linind2 = new int       [nlin + 1];
      lincoe2 = new CouNumber [nlin + 1];

      register int j;

#ifdef DEBUG
      for (j=0; j<mid;  j++) printf ("{%g x%d} ", lincoe [j], linind [j]);
      for (j++; j<nlin; j++) printf ("{%g x%d} ", lincoe [j], linind [j]);
#endif

      CouNumber divider = -1. / auxcoe;

      for (j=0; j<mid;  j++) {linind2 [j]   = linind [j]; lincoe2 [j]   = divider * lincoe [j];}
      for (j++; j<nlin; j++) {linind2 [j-1] = linind [j]; lincoe2 [j-1] = divider * lincoe [j];}

      linind2 [j-1] = -1; // terminate list of indices

#ifdef DEBUG
      for (j=0; j<mid;  j++) printf ("<%g x%d> ", lincoe2 [j],   linind2 [j]);
      for (j++; j<nlin; j++) printf ("<%g x%d> ", lincoe2 [j-1], linind2 [j-1]);
#endif

      // nl arglist is done, later decide whether to incorporate it as
      // it is or with a coefficient
    }

    // the extracted index is one term of...
    if (which >= 0) --nargs; // ...the nonlinear sum
    else            --nlin;  // ...the linear part

#ifdef DEBUG
    printf ("\n::: rhs %g, lin %d, nl %d\n", rhs, nlin, nargs);
#endif

    if ((code == COU_EXPRGROUP) && (nlin > 0)) { 

      // an exprGroup with at least one linear term left

      // build new vectors for index and coeff
      if    (fabs (auxcoe + 1) < COUENNE_EPS)

	//      f(x) + c0 -  w = rhs   =====>   w =       f(x) + c0 - rhs
	rest = new exprGroup (c0 - rhs, linind2, lincoe2, newarglist, nargs);

      else { // f(x) + c0 + aw = rhs   =====>   w = -1/a (f(x) + c0 - rhs), a != -1

	expression **mullist = new expression * [1];

	if ((nargs <= 1) && ((*newarglist) -> Linearity () <= CONSTANT)) {
	  // the only nonlinear term is a constant
	  *mullist = new exprConst ((*newarglist) -> Value ());
	  //delete *newarglist;
	  delete [] newarglist;
	}
	else // multiple nonlinear terms, multiply them by -1/auxcoe
	  *mullist = new exprMul (new exprConst (-1. / auxcoe), 
				  new exprSum (newarglist, nargs));

	// final outcome: -1/a (f(x) + c0 - rhs)
	rest = new exprGroup (-1. / auxcoe * (c0-rhs), linind2, lincoe2, mullist, 1);
      }
    }
    else { // simple exprSum

      if (fabs (rhs) > COUENNE_EPS) // have to add constant to exprSum
	if ((nargs == 1) && ((*newarglist) -> Type () == CONST)) {
	  CouNumber val = (*newarglist) -> Value () - rhs;
	  delete *newarglist;
	  *newarglist = new exprConst (val);
	} else newarglist [nargs++] = new exprConst (-rhs);

      // now exprSum is complete with -rhs. Send it to right hand side

      expression *auxDef;

      if (nargs==1) {
	auxDef = *newarglist;
	delete [] newarglist;
      } else auxDef = new exprSum (newarglist, nargs);

      if (fabs (auxcoe + 1) < COUENNE_EPS)
	rest = auxDef;
      else if ((fabs (auxcoe - 1) < COUENNE_EPS) && 
	       (auxDef -> code () == COU_EXPROPP))
	rest = auxDef -> Argument ();
      else rest = new exprMul (new exprConst (-1/auxcoe), auxDef);
    }

#ifdef DEBUG
    printf ("gotten ");
    rest -> print (); printf ("\n");
#endif

    auxInd = maxindex;

  } break;

  default: break;
  } // end switch () ////////////////////////////////////////////////////////

  if ((auxInd < 0) || (wentAux [auxInd]))
    return -1;

  // we have a variable, meaning the constraint body is of the form
  // w +/- f(x) or f(x) +/- w

  // standardize remaining of the expression

#ifdef DEBUG
  printf ("standardize rest (2nd level) "); fflush (stdout);
  rest -> print ();
#endif

  int rtype = rest -> Type ();

  if (rtype == UNARY) { //////////////////////////////////////////

    exprAux *aux = rest -> Argument () -> standardize (p);

    if (aux) {
      //delete rest -> Argument ();
      *(rest -> ArgPtr ()) = new exprClone (aux);
    }

  } else if (rtype == N_ARY) ///////////////////////////////////////////

    for (int nargs = rest -> nArgs (), i=0; i < nargs; i++) {

      exprAux *aux = rest -> ArgList () [i] -> standardize (p);

      if (aux) {
	//delete rest -> ArgList () [i];
	rest -> ArgList () [i] = new exprClone (aux);
      }
    }

  //  exprAux *aux = rest -> standardize (p);

#if 0
  printf ("... done, auxind %d\n", auxInd);

  if (aux) {

//    rest = aux -> Image () -> clone ();

//    aux -> decreaseMult ();
    printf ("... fictitious "); fflush (stdout); aux -> print (); 
    printf (" := "); aux -> Image () -> print (); 

    // find aux in vector and delete it (otherwise LP relaxation will
    // count those variables)
    for (std::vector <exprVar *>:: iterator i = p -> Variables ().begin (); 
	 i != p -> Variables ().end(); i++)

      if ((*i) -> Index() == aux -> Index ()) {
	printf (" fictitious %d to be erased: ", p -> Variables ().size ());
	(*i)             -> print (); printf (" := ");
	(*i) -> Image () -> print (); printf ("\n");
	//	p -> Variables () . erase (i);
	break;
      }
  }
#endif
  
  //  rest = aux -> Image ();

#ifdef DEBUG
  printf (" and "); fflush (stdout);
  rest -> print (); 
  //printf (", "); fflush (stdout); body -> print (); 
  printf ("\n");
#endif

  return auxInd;
}

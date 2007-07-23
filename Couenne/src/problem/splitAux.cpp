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

  case COU_EXPRSUB: { // simplest case, we have w-f(x)=rhs or f(x)-w=rhs or f(x)-g(x)=rhs

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

      pos = 1;
    }

    /*
    if ((auxInd < 0) ||                     // first argument is not a variable
	(wentAux [auxInd]))                 // or it was taken as a definition
      auxInd = alist [pos = 1] -> Index (); // default to second argument

    if ((auxInd < 0) ||     // if second argument is not a variable either
	(wentAux [auxInd])) // if it was not taken as a definition already
      return -1;            // then this is not an aux definition
    */

    //    printf ("{ind %d #%d} ", auxInd, pos); 

    expression *clone = alist [1 - pos] -> clone (); // what remains is the "independent" expression

    expression *auxdef = 
      (fabs (coeff - 1) < COUENNE_EPS) ?          // if coefficient is 1
      (clone) :                                   //     do nothing
      new exprMul (new exprConst (coeff), clone); //     else multiply it by coefficient

    rest = (fabs (rhs) < COUENNE_EPS) ?                              // no extra constant?
      (auxdef) :                                                     // just put other argument
      (new exprSum (auxdef, new exprConst ((pos==1) ? -rhs : rhs))); // otherwise sum it with \pm rhs

  } break;

    // another case is an exprGroup or exprSum, which requires to
    // decompose the expression and re-assemble

    //  case COU_EXPRQUAD: // TODO
  case COU_EXPRGROUP:
  case COU_EXPRSUM: { // more complicated, have to look for right variable

    // data structure to be used below if there is a linear term.
    // which specifies position within arrays (negative for linear
    // part of exprGroup, positive for all elements of exprSum)

    int maxindex = -1, *linind = NULL, nlin = 0, which;
    CouNumber c0 = 0., *lincoe = NULL, auxcoe = 1;

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

	  which    = - nlin - 1;    // mark which with negative number 
	  auxcoe   = lincoe [nlin];
	  maxindex = j;
	}

      //      printf (" [ind %d, coe %.1g, whi %d] ", 
      //	      maxindex, auxcoe, which);

    } //else {printf ("sum "); fflush (stdout);}

    // check indices of nonlinear elements of the sum

    //    printf (" # NL (%d) # ", body -> nArgs ());

    for (int i = body -> nArgs (); i--;) {

      //      printf ("arg %d: ", i); 

      CouNumber coeff;
      int index;

      elementBreak (alist [i], index, coeff);

      //arg -> print (); 

      if ((index > maxindex) &&
	  !(wentAux [index]) &&
	  (fabs (coeff) > COUENNE_EPS) && 
	  (body -> dependsOn (&index, 1) <= 1)) {

	maxindex = index;
	which    = i;
	auxcoe   = coeff;
      }
    }

    // if found no possible substitute, break
    if (maxindex < 0) break;

    // create a new exprGroup or exprSum with all elements but the
    // extracted auxiliary

    // start with exprSum

    int nargs = body -> nArgs ();
    expression **newarglist;

    //    printf (" [[ind %d, coe %.1g, wh %d, nargs %d, nlin %d]] ", 
    //	    maxindex, auxcoe, which, nargs, nlin);

    if (nargs > 0) { // there is an element in the nonlinear sum to be drawn

      int j, mid = (which < 0) ? nargs : which;

      newarglist = new expression * [nargs + 1];

      for (j=0; j<mid;   j++) {newarglist [j]   = alist [j]; newarglist [j]->  print ();printf(",");}
      for (j++; j<nargs; j++) {newarglist [j-1] = alist [j]; newarglist [j-1]->print ();printf(";");}

      // nl arglist is done, later decide whether to incorporate it as
      // it is or with a coefficient
    } else {

      // no nonlinear arguments, or the only one is the new aux
      newarglist  = new expression *;
      *newarglist = new exprConst (0);
    }

    auxcoe = -1. / auxcoe;

    int       *linind2 = NULL;
    CouNumber *lincoe2 = NULL;

    if (nlin > 0) { // there is an element in the linear sum to be drawn

      int mid = (which >= 0) ? nlin : - which - 1;

      linind2 = new int       [nlin + 1];
      lincoe2 = new CouNumber [nlin + 1];

      register int j;

      for (j=0; j<mid;  j++) printf ("{%g x%d} ", lincoe [j], linind [j]);
      for (j++; j<nlin; j++) printf ("{%g x%d} ", lincoe [j], linind [j]);

      for (j=0; j<mid;  j++) {linind2 [j]   = linind [j]; lincoe2 [j]   = auxcoe * lincoe [j];}
      for (j++; j<nlin; j++) {linind2 [j-1] = linind [j]; lincoe2 [j-1] = auxcoe * lincoe [j];}
      linind2 [j-1] = -1; // terminate list of indices

      for (j=0; j<mid;  j++) printf ("<%g x%d> ", lincoe2 [j],   linind2 [j]);
      for (j++; j<nlin; j++) printf ("<%g x%d> ", lincoe2 [j-1], linind2 [j-1]);

      // nl arglist is done, later decide whether to incorporate it as
      // it is or with a coefficient
    }

    // add constant(s)

    if (which >= 0) --nargs;
    else            --nlin;

    rhs = auxcoe * (rhs - c0);

    //    printf ("\n::: rhs %g, lin %d, nl %d\n", rhs, nlin, nargs);

    if ((code == COU_EXPRGROUP) && (nlin > 0)) {

      //      printf ("new group\n");
      // build new vectors for index and coeff;
      if    (fabs (auxcoe - 1) < COUENNE_EPS) 
	rest = new exprGroup (rhs, linind2, lincoe2, newarglist, nargs);
      else {

	expression **mullist = new expression * [1];

	if ((nargs == 1) && ((*newarglist) -> Linearity () <= ZERO)) {
	  *mullist = *newarglist;
	  delete [] newarglist;
	}
	else
	  *mullist = new exprMul (new exprConst (auxcoe), 
				  new exprSum (newarglist, nargs));

	rest = new exprGroup (rhs, linind2, lincoe2, mullist, 1);
      }
    }
    else {

      printf ("new sum\n");

      if (fabs (rhs) > COUENNE_EPS)
	newarglist [nargs++] = new exprConst (rhs);

      expression *auxDef;

      if (nargs==1) {
	auxDef = *newarglist;
	delete [] newarglist;
      } else auxDef = new exprSum (newarglist, nargs);

      printf ("auxdef = "); auxDef -> print ();
      printf ("\n");

      if (fabs (auxcoe - 1) < COUENNE_EPS)
	rest = auxDef;
      else if ((fabs (auxcoe + 1) < COUENNE_EPS) && 
	       (auxDef -> code () == COU_EXPROPP))
	rest = auxDef -> Argument ();
      else rest = new exprMul (new exprConst (auxcoe), auxDef);
    }

    //    printf ("gotten ");
    //    rest -> print (); printf ("\n");

    auxInd = maxindex;

    //rest = new exprGroup 
  } break;

  default: break;
  }

  if ((auxInd < 0) || (wentAux [auxInd]))
    return -1;

  // we have a variable, meaning the constraint body is of the form
  // $w \pm f(x)$ or $f(x) \pm w$

  // standardize remaining of the expression

  printf ("standardize rest "); fflush (stdout);
  rest -> print ();

  exprAux *aux = rest -> standardize (p);

  printf ("... done, auxind %d\n", auxInd);

  if (aux) {

    rest = aux -> Image () -> clone ();

    aux -> decreaseMult ();
    printf ("... fictitious "); fflush (stdout); aux -> print (); 
    printf (" := "); aux -> Image () -> print (); 

    // find aux in vector and delete it (otherwise LP relaxation will
    // count those variables)

    for (std::vector <exprVar *>:: iterator i = p -> Variables ().end(); 
	 i-- != p -> Variables ().begin();)
      if ((*i) -> Index() == aux -> Index ()) {
	printf (" fictitious %d to be erased: ", p -> Variables ().size ());
	(*i) -> print ();
	printf (" := ");
	(*i) -> Image () -> print ();
	printf ("\n");
	p -> Variables () . erase (i);
	break;
      }
  }

  //  rest = aux -> Image ();

  printf (" and "); fflush (stdout);
  rest -> print (); printf (", "); fflush (stdout);
  body -> print (); printf ("\n");

  return auxInd;
}

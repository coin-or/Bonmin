/*
 * Name:    sumStandardize.cpp
 * Author:  Pietro Belotti
 * Purpose: check if exprGroup/exprSum contains a lot of quadratic/bilinear terms
 *
 * This file is licensed under the Common Public License (CPL)
 */


#include <CouenneTypes.h>
#include <CouenneProblem.hpp>
#include <exprSum.hpp>
#include <exprMul.hpp>
#include <exprPow.hpp>
#include <exprQuad.hpp>


/// re-organizes multiplication and stores indices (and exponents) of its variables
void flattenMul (expression *mul, CouNumber &, 
		 std::map <int, CouNumber> &, 
		 CouenneProblem *);


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
int decomposeTerm (CouenneProblem *, 
		   expression *, 
		   CouNumber &, 
		   int &, int &);


/// translate a sum into:
/// 1) an exprQuad,  if some quadratic/bilinear terms exist
/// 2) an exprGroup, if only linear terms are present

exprAux *exprSum::standardize (CouenneProblem *p) {

  /*printf ("standardizing: ------------------------------\n");
  print (); 
  printf ("\n-------------------------------------------\n");*/

  // turn all elements of arglist_ and of the linear part into an exprQuad.
  // count all potential quadratic terms for exprQuad

  std::map <std::pair <int,int>, CouNumber> qmap;
  std::map <int,                 CouNumber> lmap;

  int cod = code ();

  CouNumber c0 = 0;   // final constant term

  ////////////////////////////////////////////////////////////////////////////////

  // initialize linear/quad terms with the original values/indices

  if ((cod == COU_EXPRGROUP) || 
      (cod == COU_EXPRQUAD)) {

    exprGroup *eg = dynamic_cast <exprGroup *> (this);

    c0 += eg -> getc0 ();

    int       *olind = eg -> getIndices ();
    CouNumber *olcoe = eg -> getCoeffs ();

    for (int i = eg -> getnLTerms (); i--;)
      lmap.insert (std::pair <int, CouNumber> (olind [i], olcoe [i]));

    if (cod == COU_EXPRQUAD) {

      exprQuad *eq = dynamic_cast <exprQuad *> (this);

      int       *oqindI = eq -> getQIndexI ();
      int       *oqindJ = eq -> getQIndexJ ();
      CouNumber *oqcoe  = eq -> getQCoeffs ();

      for (int i = eq -> getnQTerms (); i--;) {
	std::pair <int, int> ind (oqindI [i], oqindJ [i]);
	qmap.insert (std::pair <std::pair <int, int>, CouNumber> (ind, oqcoe [i]));
      }
    }
  }

  if (0) {
  printf ("===============================================\n");
  int i = 0;
  for (std::map <std::pair <int, int>, CouNumber>::iterator it = qmap.begin (); 
       it != qmap.end (); it++) {
    printf ("%4d %4d %5g  |  ", 
	    it -> first.first,
	    it -> first.second,
	    it -> second);

    if (!(++i % 10)) printf ("\n");
  }

  printf ("\n...................\n");

  i = 0;

  for (std::map <int, CouNumber>::iterator it = lmap.begin (); 
       it != lmap.end (); it++) {

    printf ("%4d %5g  |  ", 
	    it -> first,
	    it -> second);

    if (!(++i % 10)) printf ("\n");
  }
  printf ("\n===============================================\n");
  }

  ////////////////////////////////////////////////////////////////////////////////

  // standardize all nonlinear arguments and put them into linear or
  // quadratic form

  for (int i=0; i<nargs_; i++) {

    CouNumber coeff = 1.;
    int ind0 = -1, 
        ind1 = -1;

    //printf ("now decompose "); 
    //arglist_ [i] -> print ();

    int type = decomposeTerm (p, arglist_ [i], coeff, ind0, ind1);

    if (fabs (coeff) < COUENNE_EPS) 
      continue;

    switch (type) {

    case COU_EXPRCONST: 
      c0 += coeff; 
      break;

    case COU_EXPRVAR: {
      std::pair <int, CouNumber> npair (ind0, coeff);
      std::map <int, CouNumber>::iterator i = lmap.find (ind0);
      if (i != lmap.end()) {
	if (fabs (i -> second += coeff) < COUENNE_EPS)
	  lmap.erase (i);
      } else lmap.insert (npair);

      //lmap.insert (std::pair <int, CouNumber> (ind0, coeff));
    } break;

    case COU_EXPRMUL: {
      int i0 = ind0, i1 = ind1;
      if (i0 < i1) {i0 = ind1; i1 = ind0;}

      std::pair <int, int> nind (i0, i1);
      std::pair <std::pair <int, int>, CouNumber> npair (nind, coeff);

      std::map <std::pair <int, int>, CouNumber>::iterator i = qmap.find (nind);

      if (i != qmap.end ()) {
	if (fabs (i -> second += coeff) < COUENNE_EPS)
	  qmap.erase (i);
      } else qmap.insert (npair);

      //std::pair <int, int> ind (ind0, ind1);
      //qmap.insert (std::pair <std::pair <int, int>, CouNumber> (ind, coeff));
    } break;

    case COU_EXPRPOW: {
      std::pair <int, int> nind (ind0, ind0);
      std::pair <std::pair <int, int>, CouNumber> npair (nind, coeff);
      std::map <std::pair <int, int>, CouNumber>::iterator i = qmap.find (nind);
      if (i != qmap.end()) {
	if (fabs (i -> second += coeff) < COUENNE_EPS)
	  qmap.erase (i);
      } else qmap.insert (npair);

      //      std::pair <int, int> ind (ind0, ind0);
      //      qmap.insert (std::pair <std::pair <int, int>, CouNumber> (ind, coeff));
    } break;

    default: 
      printf ("Couenne::decomposeTerm failed\n"); 
      exit (-1);
    }
    //    printf (" --> %g,%d,%d\n", coeff, ind0, ind1);
  }

  if (0) {
  printf ("AFTER READING...\n===============================================\n");
  int i = 0;
  for (std::map <std::pair <int, int>, CouNumber>::iterator it = qmap.begin (); 
       it != qmap.end (); it++) {
    printf ("%4d %4d %5g  |  ", 
	    it -> first.first,
	    it -> first.second,
	    it -> second);

    if (!(++i % 10)) printf ("\n");
  }

  printf ("\n...................\n");

  i = 0;

  for (std::map <int, CouNumber>::iterator it = lmap.begin (); 
       it != lmap.end (); it++) {

    printf ("%4d %5g  |  ", 
	    it -> first,
	    it -> second);

    if (!(++i % 10)) printf ("\n");
  }
  printf ("\n===============================================\n");
  }

  /// data for exprQuad

  int  nq = qmap.size (),
      *qi = new int [nq+1], 
      *qj = new int [nq+1];

  CouNumber *qc = new CouNumber [nq];

  /// data for exprGroup

  int  nl = lmap.size(),
      *li = new int [nl+1];

  CouNumber *lc = new CouNumber [nl];

  // terminate arrays with negative index
  qi [nq] = li [nl] = -1; 

  std::map <int, CouNumber>::iterator lit = lmap.begin (); 

  // fill in arrays for linear part
  for (int i=0; i<nl; i++, lit++) {

    li [i] = lit -> first;
    lc [i] = lit -> second;
  }

  std::map <std::pair <int, int>, CouNumber>::iterator qit = qmap.begin (); 

  // fill in arrays for quadratic part
  for (int i=0; i < nq; i++, qit++) {
    qi [i] = qit -> first. first;
    qj [i] = qit -> first. second;
    qc [i] = qit -> second;
  }

  nl = lmap.size ();
  nq = qmap.size ();

  /*for (int i=0; i < nq;) {
    printf ("%3d %3d %4g   |  ", 
	    qi [i], qj [i], qc [i]);
    if (!(++i % 10)) printf ("\n");
  }

  printf ("\n");

  for (int i=0; i < nl;) {
    printf ("%3d %4g   |  ", 
	    li [i], lc [i]);
    if (!(++i % 10)) printf ("\n");
    }*/

  // particular cases ///////////////////////////////////////////////////////////

  // a constant
  if ((nq==0) && (nl==0))
    return p -> addAuxiliary (new exprConst (c0));

  // a linear monomial, cx
  if ((nq==0) && (fabs (c0) < COUENNE_EPS) && (nl==1)) {
    if (fabs (*lc - 1) < COUENNE_EPS) return  p -> addAuxiliary (new exprClone (p -> Var (*li)));
    else return p -> addAuxiliary (new exprMul (new exprConst (*lc), new exprClone (p -> Var (*li))));
  }

  // a bilinear/quadratic monomial, cx^2 or cxy
  if ((nl==0) && (fabs (c0) < COUENNE_EPS) && (nq==1)) {

    expression *quad;

    if (*qi == *qj) quad = new exprPow (new exprClone (p -> Var (*qi)), new exprConst (2));
    else            quad = new exprMul (new exprClone (p -> Var (*qi)), 
					new exprClone (p -> Var (*qj)));

    if (fabs (*qc - 1) < COUENNE_EPS) return p -> addAuxiliary (quad);
    else return p -> addAuxiliary (new exprMul (new exprConst (*qc), quad));
  }

  // general case ///////////////////////////////////////////////////////////////

  expression **zero = new expression * [1];
  *zero = new exprConst (0.);

  return 
    (nq==0) ? 
    (p -> addAuxiliary (new exprGroup (c0, li, lc,             zero, 1))) :
    (p -> addAuxiliary (new exprQuad  (c0, li, lc, qi, qj, qc, zero, 1)));
}

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
#include <exprSub.hpp>
#include <exprOpp.hpp>
#include <exprMul.hpp>
#include <exprPow.hpp>
#include <exprGroup.hpp>
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

void decomposeTerm (CouenneProblem *p, expression *term,
		    CouNumber initCoe,
		    CouNumber &c0,
		    std::map <int,                 CouNumber> &lmap,
		    std::map <std::pair <int,int>, CouNumber> &qmap);


exprAux *linStandardize (CouenneProblem *, CouNumber, 
			 std::map <int,                 CouNumber> &,
			 std::map <std::pair <int,int>, CouNumber> &);


/// translate a sum/difference/exprOpp into:
///
/// 1) an exprGroup, if only linear terms are present
/// 2) an exprQuad,  if some quadratic/bilinear terms exist

exprAux *exprSum::standardize (CouenneProblem *p) {

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

  ////////////////////////////////////////////////////////////////////////////////

  // standardize all nonlinear arguments and put them into linear or
  // quadratic form

  for (int i=0; i<nargs_; i++)
    decomposeTerm (p, arglist_ [i], 1, c0, lmap, qmap);

  return linStandardize (p, c0, lmap, qmap);
}



/// translate a exprOpp into:
///
/// 1) an exprGroup, if only linear terms are present
/// 2) an exprQuad,  if some quadratic/bilinear terms exist

exprAux *exprOpp::standardize (CouenneProblem *p) {

  // turn all elements of arglist_ and of the linear part into an exprQuad.
  // count all potential quadratic terms for exprQuad

  std::map <std::pair <int,int>, CouNumber> qmap;
  std::map <int,                 CouNumber> lmap;

  CouNumber c0 = 0;   // final constant term

  decomposeTerm (p, argument_, -1, c0, lmap, qmap);

  return linStandardize (p, c0, lmap, qmap);
}



/// translate a difference into:
///
/// 1) an exprGroup, if only linear terms are present
/// 2) an exprQuad,  if some quadratic/bilinear terms exist

exprAux *exprSub::standardize (CouenneProblem *p) {

  // turn all elements of arglist_ and of the linear part into an exprQuad.
  // count all potential quadratic terms for exprQuad

  std::map <std::pair <int,int>, CouNumber> qmap;
  std::map <int,                 CouNumber> lmap;

  CouNumber c0 = 0;   // final constant term

  ////////////////////////////////////////////////////////////////////////////////

  // standardize all nonlinear arguments and put them into linear or
  // quadratic form

  decomposeTerm (p, arglist_ [0],  1, c0, lmap, qmap);
  decomposeTerm (p, arglist_ [1], -1, c0, lmap, qmap);

  return linStandardize (p, c0, lmap, qmap);
}



/// standardization of linear exprOp's

exprAux *linStandardize (CouenneProblem *p, CouNumber c0, 
			 std::map <int,                 CouNumber> &lmap,
 			 std::map <std::pair <int,int>, CouNumber> &qmap) {

  ////////////////////////////////////////////////////////////////////////////////////////

  // here put sparsity analysis etc

  ////////////////////////////////////////////////////////////////////////////////////////

  int  nq = qmap.size (),     /// data for exprQuad
      *qi = new int [nq+1], 
      *qj = new int [nq+1];
  CouNumber *qc = new CouNumber [nq];

  int  nl = lmap.size(),      /// data for exprGroup
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

  exprAux *ret =
    (nq==0) ? 
    (p -> addAuxiliary (new exprGroup (c0, li, lc,             zero, 1))) :
    (p -> addAuxiliary (new exprQuad  (c0, li, lc, qi, qj, qc, zero, 1)));

  /*printf (" ==> "); 
  ret -> print (); printf (" := "); 
  ret -> Image () -> print (); printf ("\n");*/

 return ret;
}

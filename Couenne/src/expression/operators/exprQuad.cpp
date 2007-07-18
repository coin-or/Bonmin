/*
 * Name:    exprQuad.cpp
 * Author:  Pietro Belotti
 * Purpose: implementation of some methods for exprQuad
 *
 * (C) Pietro Belotti 2007. This file is licensed under the Common Public License (CPL)
 */

#include <CouenneProblem.hpp>
#include <exprConst.hpp>
#include <exprQuad.hpp>
#include <exprPow.hpp>
#include <exprMul.hpp>

/// Constructor
exprQuad::exprQuad  (CouNumber c0,      // constant term
		     int *index,        // indices (array terminated by a -1)
		     CouNumber *coeff,  // coefficient vector
		     int *qindexI,      // indices I (array terminated by a -1)
		     int *qindexJ,      // indices J (array terminated by a -1)
		     CouNumber *qcoeff, // coefficient vector
		     expression **al,   // vector of nonlinear expressions to be added 
		     int n):            // number of *nonlinear* expressions in al

  exprGroup (c0, index, coeff, al, n),
  nqterms_  (0),
  dCoeffLo_ (NULL),
  dCoeffUp_ (NULL),
  dIndex_   (NULL),
  nDiag_    (0)     {

  for (register int *qi = qindexI; *qi++ >= 0; nqterms_++);

  qindexI_ = new int       [nqterms_];
  qindexJ_ = new int       [nqterms_];
  qcoeff_  = new CouNumber [nqterms_];

  for (register int i = nqterms_; i--;) {
    qindexI_ [i] = qindexI [i];
    qindexJ_ [i] = qindexJ [i];
    qcoeff_  [i] = qcoeff  [i];
  }
} 


/// copy constructor
exprQuad::exprQuad  (const exprQuad &src): 
  exprGroup (src),

  qindexI_  (NULL),
  qindexJ_  (NULL),
  qcoeff_   (NULL),
  nqterms_  (src.nqterms_) {

  if (src.qindexI_) {
    qindexI_ = new int       [nqterms_];
    qindexJ_ = new int       [nqterms_];
    qcoeff_  = new CouNumber [nqterms_];
  }

  int *qi = src.qindexI_,
      *qj = src.qindexJ_;

  CouNumber *qc = src.qcoeff_;

  for (int i = nqterms_; i--;) {
    qindexI_ [i] = qi [i];
    qindexJ_ [i] = qj [i];
    qcoeff_  [i] = qc [i];
  }
} 


/// I/O
void exprQuad::print (std::ostream &out, bool descend, CouenneProblem *p) const {

  exprGroup::print (out, descend, p);

  for (int i = 0; i < nqterms_; i++) {

    out << qcoeff_ [i];

    if (p) {

      int qi = qindexI_ [i], 
	  qj = qindexJ_ [i];

      expression *prod;

      if (qi == qj) 
	prod    = (new exprPow (new exprClone (p -> Var (qi)), 
				new exprConst (2)));
      else prod = (new exprMul (new exprClone (p -> Var (qi)),
				new exprClone (p -> Var (qj))));

      prod -> print (out, descend, p);
      delete prod;

    } else {

      int qi = qindexI_ [i], 
	  qj = qindexJ_ [i];

      if (qi == qj) out << "x_" << qi << "^2 ";
      else          out << "x_" << qi << "*x_" << qj;
    }
  }
}


/// differentiation
expression *exprQuad::differentiate (int index) {

  expression **arglist = new expression * [nargs_+1];

  register int nonconst = 0;

  CouNumber totlin=0;
  for (register int *ind = index_, i=0; *ind>=0; i++)
    if (*ind++ == index) {
      nonconst = 1;
      totlin += coeff_ [i];
    }

  if (nonconst && (fabs (totlin) > COUENNE_EPS))
    *arglist = new exprConst (totlin);

  for (int i = 0; i < nargs_; i++) 
    if (arglist_ [i] -> dependsOn (&index, 1))
      arglist [nonconst++] = arglist_ [i] -> differentiate (index);

  if (!nonconst) {
    delete [] arglist;
    return new exprConst (0);
  }
  else return new exprSum (arglist, nonconst);
}


/// compare affine terms
int exprQuad::compare (exprQuad &e) {

  CouNumber *coe0 =   coeff_,
            *coe1 = e.coeff_;

  if (c0_ < e.c0_ - COUENNE_EPS) return -1;
  if (c0_ > e.c0_ + COUENNE_EPS) return  1;

  for (register int *ind0 = index_, *ind1 = e.index_; 
       *ind0 >= 0 || *ind1 >= 0; 
       ind0++, ind1++, 
       coe0++, coe1++) {
 
    if (*ind0 < *ind1) return -1;
    if (*ind0 > *ind1) return  1;
    if (*coe0 < *coe1 - COUENNE_EPS) return -1;
    if (*coe0 > *coe1 + COUENNE_EPS) return  1;
  }

  return 0;
}

/// used in rank-based branching variable choice

int exprQuad::rank (CouenneProblem *p) {

  int maxrank = exprOp::rank (p);

  if (maxrank < 0) 
    maxrank = 0;

  int norig = p -> nVars ();

  for (register int *ind = index_; *ind>=0; ind++) {

    int r = (*ind >= norig) ? 
      (p -> Aux (*ind - norig) -> rank (p)) :
      (p -> Var (*ind)         -> rank (p));
    if (r > maxrank)
      maxrank = r;
  }

  return maxrank;
}


/// return an index to the variable's argument that is better fixed
/// in a branching rule for solving a nonconvexity gap

expression *exprQuad::getFixVar () {
  if (arglist_ [0] -> Type () == CONST) 
    return this;
  else return arglist_ [0];
}

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

  int qi, qj;

  for (register int i = nqterms_; i--;) {
    qindexI_ [i] = qi = qindexI [i];
    qindexJ_ [i] = qj = qindexJ [i];
    qcoeff_  [i] = (qi == qj) ? (qcoeff [i]) : (0.5 * qcoeff [i]); // HIC EST DIVISION
  }
} 


/// copy constructor
exprQuad::exprQuad  (const exprQuad &src): 
  exprGroup (src),

  qindexI_  (NULL),
  qindexJ_  (NULL),
  qcoeff_   (NULL),
  nqterms_  (src.nqterms_),
  nDiag_    (src.nDiag_)   {

  if (src.qindexI_) {
    qindexI_ = new int       [nqterms_];
    qindexJ_ = new int       [nqterms_];
    qcoeff_  = new CouNumber [nqterms_];

    int *qi = src.qindexI_,
        *qj = src.qindexJ_;

    CouNumber *qc = src.qcoeff_;

    for (int i = nqterms_; i--;) {
      qindexI_ [i] = qi [i];
      qindexJ_ [i] = qj [i];
      qcoeff_  [i] = qc [i];
    }
  }

  if (src.dIndex_) {
    
    dIndex_   = new int       [nDiag_];
    dCoeffLo_ = new CouNumber [nDiag_];
    dCoeffUp_ = new CouNumber [nDiag_];

    /*
    int *qi = src.qindexI_,
        *qj = src.qindexJ_;

    CouNumber *qc = src.qcoeff_;
    */
    for (int i = nDiag_; i--;) {
      dIndex_   [i] = src.dIndex_   [i];
      dCoeffLo_ [i] = src.dCoeffLo_ [i];
      dCoeffUp_ [i] = src.dCoeffUp_ [i];
    }
  }
} 

/// I/O
void exprQuad::print (std::ostream &out, bool descend, CouenneProblem *p) const {

  // print linear and nonquadratic part
  exprGroup::print (out, descend, p);

  // print bilinear terms
  for (int i = 0; i < nqterms_; i++) {

    int qi = qindexI_ [i], 
        qj = qindexJ_ [i];

    CouNumber coe = (qi == qj) ? (qcoeff_ [i]) : (2 * qcoeff_ [i]);

    if (coe > 0) out << '+';
    out << coe;

    if (p) { // have problem pointer, use right names (x,w,y)

      expression *prod;

      if (qi == qj) 
	prod    = (new exprPow (new exprClone (p -> Var (qi)), 
				new exprConst (2)));
      else prod = (new exprMul (new exprClone (p -> Var (qi)),
				new exprClone (p -> Var (qj))));
      prod -> print (out, descend, p); out << ' ';
      delete prod;

    } else { // no problem pointer, use x for all variables

      if (qi == qj) out << "x_" << qi << "^2 ";
      else          out << "x_" << qi << "*x_" << qj << ' ';
    }
  }
}


/// differentiation
expression *exprQuad::differentiate (int index) {

 /*
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
  else return new exprSum (arglist, nonconst);*/
  // TODO!

  return NULL;
}


/// compare affine terms

int exprQuad::compare (exprQuad &e) {

  if (nqterms_ < e.nqterms_) return -1;
  if (nqterms_ > e.nqterms_) return  1;

  CouNumber *coe0 =   qcoeff_,
            *coe1 = e.qcoeff_;

  for (register int *indI0 = qindexI_, 
	            *indJ0 = qindexJ_, 
	            *indI1 = e.qindexI_, 
	            *indJ1 = e.qindexJ_, i = nqterms_; 
       i--; indI0++, indI0++, indI0++, indI0++) {
 
    if (*indI0 < *indI1) return -1;
    if (*indI0 > *indI1) return  1;

    if (*indJ0 < *indJ1) return -1;
    if (*indJ0 > *indJ1) return  1;

    if (*coe0 < *coe1 - COUENNE_EPS) return -1;
    if (*coe0 > *coe1 + COUENNE_EPS) return  1;
  }

  return 0;
}


/// used in rank-based branching variable choice

int exprQuad::rank (CouenneProblem *p) {

  int maxrank = exprGroup::rank (p);

  if (maxrank < 0) 
    maxrank = 0;

  CouNumber *coe = qcoeff_;

  int  n = nqterms_, 
      *i = qindexI_,
      *j = qindexJ_;

  while (n--) 
    if (fabs (*coe++) > COUENNE_EPS) {

      register int r;

      if ((r = p -> Var (*i) -> rank (p)) > maxrank) maxrank = r;
      if ((r = p -> Var (*j) -> rank (p)) > maxrank) maxrank = r;
    }

  return maxrank;
}


/// return an index to the variable's argument that is better fixed
/// in a branching rule for solving a nonconvexity gap

expression *exprQuad::getFixVar () {

  // TODO: this is quite complicated. It is a nonlinear expression but
  // we have no access to variable pointers
  if (arglist_ [0] -> Type () == CONST) 
    return this;
  else return arglist_ [0];
}

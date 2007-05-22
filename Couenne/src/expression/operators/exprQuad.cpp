/*
 * Name:    exprQuad.cpp
 * Author:  Pietro Belotti
 * Purpose: implementation of some methods for exprQuad
 *
 * (C) Pietro Belotti 2007. This file is licensed under the Common Public License (CPL)
 */

#include <CouenneProblem.h>
#include <exprConst.h>
#include <exprQuad.h>

/// Constructor
exprQuad::exprQuad  (CouNumber c0,     // constant term
		       int *index,       // indices (array terminated by a -1)
		       CouNumber *coeff, // coefficient vector
		       int *qindexI,       // indices I (array terminated by a -1)
		       int *qindexJ,       // indices J (array terminated by a -1)
		       CouNumber *qcoeff, // coefficient vector
		       expression **al,  // vector of nonlinear expressions to be added 
		       int n):           // number of *nonlinear* expressions in al
  exprGroup (c0, index, coeff, al, n) {
  //  c0_     (c0) {
  /*  
  int nlin = 0;

  // count linear terms
  for (register int *ind = index; *ind++ >= 0; nlin++);

  index_ = new int       [nlin + 1];
  coeff_ = new CouNumber [nlin];

  index_ [nlin] = index [nlin];

  while (nlin--) {
    index_ [nlin] = index [nlin];
    coeff_ [nlin] = coeff [nlin];
  }
  */
} 


/// copy constructor
exprQuad::exprQuad  (const exprQuad &src): 
  exprGroup (src),
  //  c0_     (src.c0_),
  qindexI_  (NULL),
  qindexJ_  (NULL),
  qcoeff_  (NULL)  {
  /*
  register int *ind, size = 0;

  for (ind = src.index_; (*ind++) >= 0; size++);

  // ind is now PAST the -1, put it back on it
  --ind;

  coeff_ = new CouNumber [size];
  index_ = new int       [size+1];

  index_ [size] = -1;

  while (size--) {

    index_ [size] = *--ind;
    coeff_ [size] = src.coeff_ [size];
  }
  */
} 


/// I/O
void exprQuad::print (std::ostream &out, bool descend, CouenneProblem *p) const {
//void exprQuad::print (std::ostream &out, CouenneProblem *p = NULL) const {

//  if (nargs_ && ((nargs_ > 1) ||
//		 ((*arglist_) -> Type () != CONST) ||
//		 (fabs ((*arglist_) -> Value ()) > COUENNE_EPS)))
  exprGroup::print (out, descend, p);
  /*
  int nOrig = p ? (p -> nVars ()) : -1;

  if      (c0_ >   COUENNE_EPS) out << '+' << c0_;
  else if (c0_ < - COUENNE_EPS) out        << c0_;

  for (register int *ind=index_, i=0; *ind>=0; ind++) {

    CouNumber coeff = coeff_ [i++];

    out << ' ';

    if      (coeff >   COUENNE_EPS) out << '+' << coeff << "*";
    else if (coeff < - COUENNE_EPS) out        << coeff << "*";
    else continue;

    if (nOrig < 0) out << "x_" << *ind;
    else {
      //      out << "(";
      if (*ind < nOrig) p -> Var (*ind)       -> print (out, descend, p);
      else              p -> Aux (*ind-nOrig) -> print (out, descend, p);
      //      out << ")";
    }
    }*/

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


/// get a measure of "how linear" the expression is:
int exprQuad::Linearity () {
  /*
  // compute linearity of nonlinear part
  int nllin = exprSum::Linearity ();

  if (nllin == ZERO) // if nonlinear part equals zero
    if (*index_ == -1)
      if (fabs (c0_) < COUENNE_EPS) return ZERO;
      else                          return CONSTANT; 
    else                            return LINEAR;
  else // if nonlinear part is anything but zero
    if (nllin == CONSTANT)
      if (*index_ == -1)            return CONSTANT; 
      else                          return LINEAR;
    else                            return nllin;
  */
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
    if (++r > maxrank) // increment because above exprOp::rank returns
		       // something already incremented
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

/*
 * Name:    exprGroup.cpp
 * Author:  Pietro Belotti
 * Purpose: implementation of some methods for exprGroup
 *
 * (C) Pietro Belotti 2007. This file is licensed under the Common Public License (CPL)
 */

#include <CouenneProblem.hpp>
#include <exprConst.hpp>
#include <exprGroup.hpp>
#include <depGraph.hpp>

/// Constructor
exprGroup::exprGroup  (CouNumber c0,     // constant term
		       int *index,       // indices (array terminated by a -1)
		       CouNumber *coeff, // coefficient vector
		       expression **al,  // vector of nonlinear expressions to be added 
		       int n):           // number of *nonlinear* expressions in al
  exprSum  (al, n),
  c0_      (c0),
  nlterms_ (0) {


  // count linear terms
  for (register int *ind = index; *ind++ >= 0; nlterms_++);

  index_ = new int       [nlterms_ + 1];
  coeff_ = new CouNumber [nlterms_];

  index_ [nlterms_] = index [nlterms_]; // assign -1 at end of array

  for (int i=0; i<nlterms_; i++) {
    index_ [i] = index [i];
    coeff_ [i] = coeff [i];
  }
} 


/// copy constructor
exprGroup::exprGroup  (const exprGroup &src): 
  exprSum   (src.clonearglist (), src.nargs_),
  c0_       (src.c0_),
  index_    (NULL),
  coeff_    (NULL),
  nlterms_  (src.nlterms_) {

  coeff_ = new CouNumber [nlterms_];
  index_ = new int       [nlterms_ + 1];

  index_ [nlterms_] = -1;

  for (int i=0; i < nlterms_; i++) {

    index_ [i] = src.index_ [i];
    coeff_ [i] = src.coeff_ [i];
  }
} 


/// I/O
void exprGroup::print (std::ostream &out, bool descend, CouenneProblem *p) const {

  if (nargs_ && ((nargs_ > 1) ||
		 ((*arglist_) -> Type () != CONST) ||
		 (fabs ((*arglist_) -> Value ()) > COUENNE_EPS)))
    exprSum::print (out, descend, p);

  if      (c0_ >   COUENNE_EPS) out << '+' << c0_;
  else if (c0_ < - COUENNE_EPS) out        << c0_;

  for (register int *ind=index_, i=0; *ind>=0; ind++) {

    CouNumber coeff = coeff_ [i++];
    out << ' ';

    if      (coeff >   COUENNE_EPS) out << '+' << coeff << "*";
    else if (coeff < - COUENNE_EPS) out        << coeff << "*";
    else continue;

    if (!p) out << "x_" << *ind;
    else p -> Var (*ind) -> print (out, descend, p);
  }
}


/// differentiation
expression *exprGroup::differentiate (int index) {

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
int exprGroup::Linearity () {

  // compute linearity of nonlinear part
  int nllin = exprSum::Linearity ();

  int llin  = (nlterms_ == 0) ? 
    ((fabs (c0_) < COUENNE_EPS) ? ZERO : CONSTANT) : 
    LINEAR;

  return (llin > nllin) ? llin : nllin;
}


/// compare affine terms
int exprGroup::compare (exprGroup &e) {

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

int exprGroup::rank (CouenneProblem *p) {

  int maxrank = exprOp::rank (p);

  if (maxrank < 0) 
    maxrank = 0;

  //int norig = p -> nVars ();

  for (register int *ind = index_; *ind>=0; ind++) {

    int r = //(*ind >= norig) ? 
      //(p -> Aux (*ind - norig) -> rank (p)) :
      (p -> Var (*ind)         -> rank (p));

    if (r > maxrank)
      maxrank = r;
  }

  return maxrank;
}


/// check if this expression depends on a set of variables specified
/// in the parameters
int exprGroup::dependsOn (register int *chg, register int nch) {

  int tot = exprOp::dependsOn (chg, nch);

  /// TODO: check if chg and index_ are sorted to improve efficiency
  for (; nch-- && (*chg >= 0); chg++)
    for (register int *ind = index_; *ind >= 0;)
      if (*ind++ == *chg) 
	tot++;

  return tot;
}


/// update dependence set with index of this variable
void exprGroup::fillDepSet (std::set <DepNode *, compNode> *dep, DepGraph *g) {

  exprOp::fillDepSet (dep, g);

  for (int *index = index_; *index >= 0; index++)
    dep -> insert (g -> lookup (*index));
}

/// specialized version to check expression of linear term
int exprGroup::dependsOn (CouenneProblem *p, int *chg, int nch) {

  int tot = dependsOn (chg, nch);

  for (int *ind = index_; *ind >= 0; ind++)
    for (int i=0; (i < nch) && (chg [i] >= 0); i++)
      if (p -> Var (*ind) -> Type () == AUX)
	tot += p -> Var (*ind) -> Image () -> dependsOn (p, chg, nch);

  return tot;
}


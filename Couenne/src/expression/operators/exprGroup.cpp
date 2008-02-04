/*
 * Name:    exprGroup.cpp
 * Author:  Pietro Belotti
 * Purpose: implementation of some methods for exprGroup
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "exprConst.hpp"
#include "exprVar.hpp"
#include "exprGroup.hpp"
#include "depGraph.hpp"

class Domain;

/// Constructor
exprGroup::exprGroup  (CouNumber c0,
		       std::vector <std::pair <exprVar *, CouNumber> > &lcoeff, 
		       expression **al, 
		       int n):
  exprSum  (al, n),
  lcoeff_  (lcoeff),
  c0_      (c0) {}


/// copy constructor
exprGroup::exprGroup  (const exprGroup &src, Domain *d): 
  exprSum   (src.clonearglist (d), src.nargs_),
  c0_       (src.c0_) {

  for (lincoeff::iterator i = src.lcoeff_.begin (); i != src.lcoeff_.end (); ++i)
    lcoeff_ . push_back (std::pair <exprVar *, CouNumber> 
			 (dynamic_cast <exprVar *> (i -> first -> clone (d)), i -> second));
}


/// I/O
void exprGroup::print (std::ostream &out, bool descend) const {

  if (nargs_ && ((nargs_ > 1) ||
		 ((*arglist_) -> Type () != CONST) ||
		 (fabs ((*arglist_) -> Value ()) > COUENNE_EPS)))
    exprSum::print (out, descend);

  if      (c0_ >   0.) out << '+' << c0_;
  else if (c0_ < - 0.) out        << c0_;

  for (int n = lcoeff_.size (), i=0; n--; i++) {

    CouNumber coeff = lcoeff_ [i]. second;
    out << ' ';

    if      (coeff >   0.) out << '+' << coeff << "*";
    else if (coeff < - 0.) out        << coeff << "*";
    //else continue;

    lcoeff_ [i]. first -> print (out, descend);
  }
}


/// differentiation
expression *exprGroup::differentiate (int index) {

  expression **arglist = new expression * [nargs_ + 1];

  CouNumber totlin=0;

  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el)
    if (el -> first -> Index () == index)
      totlin += el -> second;

  int nargs = 0;

  if (fabs (totlin) > COUENNE_EPS)
    arglist [nargs++] = new exprConst (totlin);

  for (int i = 0; i < nargs_; i++) 
    if (arglist_ [i] -> dependsOn (&index, 1))
      arglist [nargs++] = arglist_ [i] -> differentiate (index);

  if ((nargs == 0) ||
      (nargs == 1) && (fabs (totlin) > COUENNE_EPS)) {
    delete [] arglist;
    return new exprConst (totlin);
  }
  else return new exprSum (arglist, nargs);
}


/// get a measure of "how linear" the expression is:
int exprGroup::Linearity () {

  int 
    nllin = exprSum::Linearity (),    // linearity of nonlinear part
    llin  = (lcoeff_.size () == 0) ?  //              linear part
    ((fabs (c0_) < COUENNE_EPS) ? ZERO : CONSTANT) : 
    LINEAR;

  return (llin > nllin) ? llin : nllin;
}


/// compare affine terms
int exprGroup::compare (exprGroup &e) {

  //int sum = exprSum::compare (e);

  //if (sum != 0) 
  //return sum;

  if (c0_ < e.c0_ - COUENNE_EPS) return -1;
  if (c0_ > e.c0_ + COUENNE_EPS) return  1;

  if (lcoeff_.size () < e.lcoeff_.size ()) return -1;
  if (lcoeff_.size () > e.lcoeff_.size ()) return  1;

  for (lincoeff::iterator 
	 el1 =   lcoeff_.begin (),
	 el2 = e.lcoeff_.begin ();
       el1 != lcoeff_.end (); 
       ++el1, ++el2) {

    int 
      ind1 = el1 -> first -> Index (),
      ind2 = el2 -> first -> Index ();

    CouNumber 
      coe1 = el1 -> second,
      coe2 = el2 -> second;

    if (ind1 < ind2) return -1;
    if (ind1 > ind2) return  1;

    if (coe1 < coe2 - COUENNE_EPS) return -1;
    if (coe1 > coe2 + COUENNE_EPS) return  1;
  }

  return 0;
}


/// used in rank-based branching variable choice

int exprGroup::rank () {

  int maxrank = exprOp::rank ();

  if (maxrank < 0) 
    maxrank = 0;

  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el) {

    int r = el -> first -> rank ();
    if (r > maxrank)
      maxrank = r;
  }

  return maxrank;
}


/// update dependence set with index of this variable
void exprGroup::fillDepSet (std::set <DepNode *, compNode> *dep, DepGraph *g) {

  exprOp::fillDepSet (dep, g);

  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el)
    dep -> insert (g -> lookup (el -> first -> Index ()));
}


/// fill in the set with all indices of variables appearing in the
/// expression
int exprGroup::DepList (std::set <int> &deplist,
			enum dig_type type) {

  int deps = exprOp::DepList (deplist, type);

  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el)
    deps += el -> first -> DepList (deplist, type);

  return deps;
}


/// is this linear term integer?
bool exprGroup::isInteger () {

  if (!(::isInteger (c0_)) ||
      !(exprOp::isInteger ()))
    return false;

  for (lincoeff::iterator el = lcoeff_.begin (); el != lcoeff_.end (); ++el) {

    CouNumber coe = el -> second;

    bool
      intCoe = ::isInteger (coe),
      intVar = el -> first -> isInteger ();

    if (intCoe && intVar)
      continue;

    CouNumber lb = (*(el -> first -> Lb ())) ();

    // check var fixed and product is integer
    if ((fabs (lb - (*(el -> first -> Ub ())) ()) < COUENNE_EPS) &&
	(::isInteger (lb * coe) ||
	 (intCoe && ::isInteger (lb)))) 
      continue;

    return false;
  }

  return true;
}


/// replace variable x with new (aux) w
void exprGroup::replace (exprVar *x, exprVar *w) {

  exprOp::replace (x, w);

  int index = x -> Index ();

  for (lincoeff::iterator el = lcoeff_. begin (); el != lcoeff_. end (); ++el) {

    exprVar * &ve = el -> first;

    if ((ve -> Type  () == VAR) &&
	(ve -> Index () == index))
      ve = w;
  }
}


/// return pointer to variable domain
Domain *exprGroup::domain () {
  if (lcoeff_.size () > 0)
    return lcoeff_ [0]. first -> domain ();
  return exprOp::domain ();
}

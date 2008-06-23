/**
 * Name:    lqelems.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of elemental elements of linear and bilinear expressions
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_LQELEMS_H
#define COUENNE_LQELEMS_H

#include <map>

#include "CouenneTypes.hpp"

class quadElem {

private:
  exprVar   *varI_;
  exprVar   *varJ_;
  CouNumber  coeff_;

public:

  quadElem (exprVar *i, exprVar *j, CouNumber c):
    varI_ (i),
    varJ_ (j),
    coeff_ (c) {}

  quadElem (const quadElem &src):
    varI_ (src.varI_),
    varJ_ (src.varJ_),
    coeff_ (src.coeff_) {}

  quadElem *clone ()
  {return new quadElem (*this);}

  inline exprVar   *varI  () {return varI_;}
  inline exprVar   *varJ  () {return varJ_;}
  inline CouNumber  coeff () {return coeff_;}
};


class LinMap {

private:
  std::map <int, CouNumber> lmap_;

public:

  /// public access
  std::map <int, CouNumber> &Map () 
  {return lmap_;}

  /// insert a pair <int,CouNumber> into a map for linear terms
  void insert (int index, CouNumber coe) {

    std::map <int, CouNumber>::iterator i = lmap_.find (index);

    if (i != lmap_.end()) {
      if (fabs (i -> second += coe) < COUENNE_EPS)
	lmap_.erase (i);
    } else {
      std::pair <int, CouNumber> npair (index, coe);
      lmap_.insert (npair);
    }
  }
};


class QuadMap {

private:
  std::map <std::pair <int, int>, CouNumber> qmap_;

public:

  /// public access
  std::map <std::pair <int, int>, CouNumber> &Map () 
  {return qmap_;}

  /// insert a pair <<int,int>,CouNumber> into a map for quadratic terms
  void insert (int indI, int indJ, CouNumber coe) {

    std::pair <int, int> nind (indI, indJ);
    std::map  <std::pair <int, int>, CouNumber>::iterator i = qmap_.find (nind);

    if (i != qmap_.end()) {
      if (fabs (i -> second += coe) < COUENNE_EPS)
	qmap_.erase (i);
    } else {
      std::pair <std::pair <int, int>, CouNumber> npair (nind, coe);
      qmap_.insert (npair);
    }
  }
};

#endif

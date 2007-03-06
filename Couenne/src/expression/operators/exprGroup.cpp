/*
 * Name:    exprGroup.cpp
 * Author:  Pietro Belotti
 * Purpose: implementation of some methods for exprGroup
 *
 * (C) Pietro Belotti 2007. This file is licensed under the Common Public License (CPL)
 */

#include <exprGroup.h>

/// copy constructor
exprGroup::exprGroup  (const exprGroup &src): 
  exprSum (clonearglist (), src.nargs_),
  c0_     (src.c0_),
  index_  (NULL),
  coeff_  (NULL)  {

  register int *ind;

  for (ind = src.index_; (*ind++) >= 0;);

  register int size = ind - src.index_;
  coeff_ = new CouNumber [size];
  index_ = new int       [size];

  index_ [--size] = -1;

  while (size--) {

    index_ [size] = *ind--;
    coeff_ [size] = src.coeff_ [size];
  }
} 


/// String equivalent (for comparisons)
const std::string exprGroup::name () const {

}


/// I/O
void exprGroup::print (std::ostream &) const {

}


/// differentiation
expression *exprGroup::differentiate (int index) {

}


/// get a measure of "how linear" the expression is:
int exprGroup::Linearity () {

}

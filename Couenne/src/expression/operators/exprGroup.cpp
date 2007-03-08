/*
 * Name:    exprGroup.cpp
 * Author:  Pietro Belotti
 * Purpose: implementation of some methods for exprGroup
 *
 * (C) Pietro Belotti 2007. This file is licensed under the Common Public License (CPL)
 */

#include <exprConst.h>
#include <exprGroup.h>

/// Constructor
exprGroup::exprGroup  (CouNumber c0, 
		       int *index, CouNumber *coeff, 
		       expression **al, int n): 
  exprSum (al, n),
  c0_     (c0) {
  
  int nlin = 0;

  for (register int *ind = index; *ind++>=0; nlin++);

  index_ = new int       [nlin+1];
  coeff_ = new CouNumber [nlin];

  index_ [nlin] = index [nlin];

  while (nlin--) {
    index_ [nlin] = index [nlin];
    coeff_ [nlin] = coeff [nlin];
  }
} 


/// copy constructor
exprGroup::exprGroup  (const exprGroup &src): 
  exprSum (src.clonearglist (), src.nargs_),
  c0_     (src.c0_),
  index_  (NULL),
  coeff_  (NULL)  {

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
} 

#define MAX_NAME 10000

/// String equivalent (for comparisons)
const std::string exprGroup::name () const {

  std::string nl = exprSum::name ();

  register char *s   = (char *) malloc (MAX_NAME * sizeof (char));
  char          *sub = (char *) malloc (MAX_NAME * sizeof (char));

  sprintf (s, "%s", nl. c_str ());

  if (fabs (c0_) > COUENNE_EPS) {
    sprintf (sub, "+%f", c0_);
    strcat (s, sub);
  }

  for (register int *ind=index_, i=0; *ind>=0;) {
    sprintf (sub, "%+fx_%d", coeff_ [i++], *ind++);
    strcat (s, sub);
  }

  s = (char *) realloc (s, (1 + strlen (s)) * sizeof (char));
  std::string ret (1 + strlen (s), ' ');
  for (register int i=strlen (s); i--;)
    ret [i] = s [i]; 

  free (s);
  free (sub);

  return ret;
}


/// I/O
void exprGroup::print (std::ostream &out) const {

  exprSum::print (out);

  if      (c0_ >   COUENNE_EPS) out << '+' << c0_;
  else if (c0_ < - COUENNE_EPS) out        << c0_;

  for (register int *ind=index_, i=0; *ind>=0;) {

    CouNumber coeff = coeff_ [i++];

    out << ' ';

    if      (coeff >   COUENNE_EPS) out << '+' << coeff;
    else if (coeff < - COUENNE_EPS) out        << coeff;

    out << " x_" << *ind++;
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
}

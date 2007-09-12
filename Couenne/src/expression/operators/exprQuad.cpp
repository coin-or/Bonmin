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
#include <depGraph.hpp>


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

  // count quadratic terms
  for (register int *qi = qindexI; *qi++ >= 0; nqterms_++);

  qindexI_ = new int       [nqterms_];
  qindexJ_ = new int       [nqterms_];
  qcoeff_  = new CouNumber [nqterms_];

  int qi, qj;

  for (register int i = nqterms_; i--;) {

    qindexI_ [i] = qi = qindexI [i];
    qindexJ_ [i] = qj = qindexJ [i];
    qcoeff_  [i] = (qi == qj) ? (qcoeff [i]) : (0.5 * qcoeff [i]); // Division
  }
} 


/// copy constructor
exprQuad::exprQuad  (const exprQuad &src): 
  exprGroup (src),

  qindexI_  (NULL),
  qindexJ_  (NULL),
  qcoeff_   (NULL),
  nqterms_  (src.nqterms_),
  dCoeffLo_ (NULL),
  dCoeffUp_ (NULL),
  dIndex_   (NULL), 
  nDiag_    (src.nDiag_) {

  // copy quadratic part (if any)

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

  // copy convexification part (if any)
  
  if (src.dIndex_) {
    
    dIndex_                     = new int       [nDiag_];
    dCoeffLo_ = (src.dCoeffLo_) ? new CouNumber [nDiag_] : NULL; 
    dCoeffUp_ = (src.dCoeffUp_) ? new CouNumber [nDiag_] : NULL; 

    for                    (int i = nDiag_; i--;) dIndex_   [i] = src.dIndex_   [i];
    if (src.dCoeffLo_) for (int i = nDiag_; i--;) dCoeffLo_ [i] = src.dCoeffLo_ [i];
    if (src.dCoeffUp_) for (int i = nDiag_; i--;) dCoeffUp_ [i] = src.dCoeffUp_ [i];
  }
} 


/// Destructor
exprQuad::~exprQuad () {

  register expression *elem;

  if (arglist_) {
    for (register int i = nargs_; i--;)
      if ((elem = arglist_ [i]))
	delete elem;

    delete [] arglist_;
    arglist_ = NULL;
  }

  if (qindexI_) {
    delete [] qindexI_;
    delete [] qindexJ_;
    delete [] qcoeff_;
  }

  if (dIndex_) {
    delete [] dIndex_;
    if (dCoeffLo_) delete [] dCoeffLo_;
    if (dCoeffUp_) delete [] dCoeffUp_;
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
    out << coe << '*';

    if (p) { // have problem pointer, use right names (x,w,y)

      expression *prod;

      if (qi == qj) 
	prod    = (new exprPow (new exprClone (p -> Var (qi)), 
				new exprConst (2)));
      else prod = (new exprMul (new exprClone (p -> Var (qi)),
				new exprClone (p -> Var (qj))));
      prod -> print (out, descend, p); out << ' ';
      delete prod;

    } else { // there is no problem pointer, use x for all variables

      if (qi == qj) out << "x_" << qi << "^2 ";
      else          out << "x_" << qi << "*x_" << qj << ' ';
    }
  }
}


/// differentiation
expression *exprQuad::differentiate (int index) {

  std::map <int, CouNumber> lmap;

  CouNumber c0 = 0;

  // derive linear part (obtain constant)
  for (register int *ind = index_, i=0; *ind>=0; i++)
    if (*ind++ == index)
      c0 += coeff_ [i];

  // derive quadratic part (obtain linear part)
  for (register int *qi = qindexI_, *qj = qindexJ_, i=0; 
       i < nqterms_; i++, qi++, qj++)

    if      (*qi == index)
      if    (*qj == index) linsert (lmap, index, 2 * qcoeff_ [i]);
      else                 linsert (lmap, *qj,       qcoeff_ [i]);
    else if (*qj == index) linsert (lmap, *qi,       qcoeff_ [i]);

  // derive nonlinear sum
  expression **arglist = new expression * [nargs_ + 1];
  int nargs = 0;

  for (int i = 0; i < nargs_; i++) 
    if (arglist_ [i] -> dependsOn (&index, 1))
      arglist [nargs++] = arglist_ [i] -> differentiate (index);

  // special cases

  // 1) no linear part
  if (lmap.empty ()) {

    // and no nonlinear part either
    if (!nargs) {
      delete arglist;
      return new exprConst (c0);
    }

    if (fabs (c0) > COUENNE_EPS)
      arglist [nargs++] = new exprConst (c0);

    return new exprSum (arglist, nargs);
  }

  // translate lmap into a vector

  int nl = lmap.size(), *linin = new int [1 + nl], j = 0;
  CouNumber *coeff = new CouNumber [nl];

  for (std::map <int, CouNumber>::iterator i = lmap.begin (); i != lmap.end (); i++) {
    linin [j]   = i -> first;
    coeff [j++] = i -> second;
  }

  linin [j] = -1;

  return new exprGroup (c0, linin, coeff, arglist, nargs);
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

  return NULL;

  // TODO: this is quite complicated. It is a nonlinear expression but
  // we have no access to variable pointers

  //if (arglist_ [0] -> Type () == CONST) 
  //  return this;
  //else return arglist_ [0];
}


/// update dependence set with index of this variable
void exprQuad::fillDepSet (std::set <DepNode *, compNode> *dep, DepGraph *g) {

  exprGroup::fillDepSet (dep, g);

  for (int *qi = qindexI_, *qj = qindexJ_, n = nqterms_; n--;) {
    dep -> insert (g -> lookup (*qi++));
    dep -> insert (g -> lookup (*qj++));
  }
}


/// insert a pair <int,CouNumber> into a map for linear terms
void linsert (std::map <int, CouNumber> &lmap, 
	      int index, CouNumber coe) {

  std::map <int, CouNumber>::iterator i = lmap.find (index);
  if (i != lmap.end()) {
    if (fabs (i -> second += coe) < COUENNE_EPS)
      lmap.erase (i);
  } else {
    std::pair <int, CouNumber> npair (index, coe);
    lmap.insert (npair);
  }
}

/// insert a pair <<int,int>,CouNumber> into a map for quadratic terms
void qinsert (std::map <std::pair <int, int>, CouNumber> &map, 
	      int indI, int indJ, CouNumber coe) {

  std::pair <int, int> nind (indI, indJ);
  std::pair <std::pair <int, int>, CouNumber> npair (nind, coe);
  std::map  <std::pair <int, int>, CouNumber>::iterator i = map.find (nind);

  if (i != map.end()) {
    if (fabs (i -> second += coe) < COUENNE_EPS)
      map.erase (i);
  } else map.insert (npair);
}


/// create dIndex_ based on occurrences in qindexI_ and qindexJ_
void exprQuad::make_dIndex (int numcols, int *indexmap) { 

  int *qindexI = getQIndexI (),
      *qindexJ = getQIndexJ ();

  nDiag_ = 0; // initialize sparse diagonal index vector

  // fill indexmap
  for (int i = 0; i < nqterms_; ++i) {

    int qi = qindexI [i], 
        qj = qindexJ [i];

    if                (indexmap [qi] == -1)  indexmap [qi] = nDiag_++;
    if ((qi != qj) && (indexmap [qj] == -1)) indexmap [qj] = nDiag_++;
  }

  dIndex_ = new int [nDiag_];

  for (int i=0; i<numcols; ++i)
    if (indexmap [i] > -1)
      dIndex_ [indexmap [i]] = i;
}

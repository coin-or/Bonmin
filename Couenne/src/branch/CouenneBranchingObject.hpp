/*
 * Name:    CouenneBranchingObject.hpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Branching object for auxiliary variables
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNEBRANCHINGOBJECT_HPP
#define COUENNEBRANCHINGOBJECT_HPP

#include <CoinFinite.hpp>
#include <OsiBranchingObject.hpp>
#include <exprAux.hpp>

#define COUENNE_CROP 1
#define COUENNE_LCROP (1e2*COUENNE_CROP)

#define COUENNE_LARGE_INTERVAL 1e4
#define COUENNE_NEAR_BOUND 1e-2


/* "Spatial" branching object. 
 *
 *  Branching can also be performed on continuous variables.
 */

class CouenneBranchingObject: public OsiTwoWayBranchingObject {

public:

  /// Return global value for convex combination between current point
  /// and midpoint.
  static CouNumber Alpha () 
  {return alpha_;}

  /// Constructor
  CouenneBranchingObject (int, int, CouNumber = - COIN_DBL_MAX, bool = false);

  /// Copy constructor
  CouenneBranchingObject (const CouenneBranchingObject &src):
    OsiTwoWayBranchingObject (src),
    index_   (src.index_),
    integer_ (src.integer_) {}

  /// Cloning method
  virtual OsiBranchingObject * clone () const
  {return new CouenneBranchingObject (*this);}

  /** \brief Execute the actions required to branch, as specified by the
	     current state of the branching object, and advance the object's
	     state. 
	     Returns change in guessed objective on next branch
  */
  virtual double branch (OsiSolverInterface * solver = NULL);

protected:

  /// The index of the variable this branching object refers to. If
  /// the corresponding CouenneObject was created on w=f(x,y), it is
  /// either x or y, chosen previously with a call to getFixVar()
  /// expression *reference_;
  int index_;

  /// True if the related variable is integer
  bool integer_;

  /// Global value for convex combination between current point and
  /// midpoint
  static CouNumber alpha_;
};

/// returns a point "inside enough" a given interval, or x if it is already
CouNumber midInterval (CouNumber x, CouNumber l, CouNumber u);

#endif

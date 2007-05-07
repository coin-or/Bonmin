/*
 * Name:    CouenneThreeWayBranchObj.hpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Three way branching object for auxiliary variables
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNETHREEWAYBRANCHOBJ_HPP
#define COUENNETHREEWAYBRANCHOBJ_HPP

#include <CoinFinite.hpp>
#include <OsiBranchingObject.hpp>
#include <exprAux.h>


/// Spatial, three-way branching object. Branching is performed on
/// continuous variables but a better convexification is sought around
/// the current point by dividing the interval in three parts

class CouenneThreeWayBranchObj: public OsiBranchingObject {

public:

  /// Constructor
  CouenneThreeWayBranchObj (expression *, CouNumber, CouNumber, CouNumber);

  /// Copy constructor
  CouenneThreeWayBranchObj (const CouenneThreeWayBranchObj &src):
    OsiBranchingObject (src),
    reference_ (src.reference_),
    lcrop_     (src.lcrop_),
    rcrop_     (src.rcrop_) {}

  /// Cloning method
  virtual OsiBranchingObject * clone () const
  {return new CouenneThreeWayBranchObj (*this);}

  /** \brief Execute the actions required to branch, as specified by the
	     current state of the branching object, and advance the object's
	     state. 
	     Returns change in guessed objective on next (what does
	     "next" mean here?) branch
  */
  virtual double branch (OsiSolverInterface * solver = NULL);

protected:

  /// The variable this branching object refers to. If the
  /// corresponding CouenneObject was created on w=f(x,y), it is
  /// either x or y.
  expression *reference_;

  /// Dividing points of interval
  CouNumber lcrop_;
  CouNumber rcrop_;

  /// first branch to be performed: 0 is left, 1 is central, 2 is right
  int first_;
};

#endif

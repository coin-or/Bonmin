/*
 * Name:    CouenneBranchingObject.hpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Branching object for auxiliary variables
 *
 * (C) Carnegie-Mellon University, 2006-07. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNEBRANCHINGOBJECT_HPP
#define COUENNEBRANCHINGOBJECT_HPP

#include "CoinFinite.hpp"
#include "OsiBranchingObject.hpp"
#include "exprAux.hpp"
#include "CouenneJournalist.hpp"

#define COUENNE_CROP 1
#define COUENNE_LCROP (1e2*COUENNE_CROP)

#define COUENNE_LARGE_INTERVAL 1e4
#define COUENNE_NEAR_BOUND 1e-2

// use to test how branching point moves. Define as index of
// independent variable in log expression
#define BR_TEST_LOG -1
#define BR_TEST_GRAPH 0

/** "Spatial" branching object. 
 *
 *  Branching can also be performed on continuous variables.
 */

class CouenneBranchingObject: public OsiTwoWayBranchingObject {

public:

  /// Constructor
  CouenneBranchingObject (JnlstPtr jnlst, expression *, int, CouNumber = - COIN_DBL_MAX);

  /// Copy constructor
  CouenneBranchingObject (const CouenneBranchingObject &src):
    OsiTwoWayBranchingObject (src),
    variable_     (src.variable_),
    jnlst_        (src.jnlst_),
    doFBBT_       (src.doFBBT_),
    doConvCuts_   (src.doConvCuts_),
    downEstimate_ (src.downEstimate_),
    upEstimate_   (src.upEstimate_) {}

  /// Cloning method
  virtual OsiBranchingObject * clone () const
  {return new CouenneBranchingObject (*this);}

  /** \brief Execute the actions required to branch, as specified by the
	     current state of the branching object, and advance the object's
	     state. 
	     Returns change in guessed objective on next branch
  */
  virtual double branch (OsiSolverInterface * solver = NULL);

  /// does this branching object only change variable bounds?
  virtual bool boundBranch () const
  {return !doConvCuts_;} // no, if it adds convexification cuts

protected:

  /// The index of the variable this branching object refers to. If
  /// the corresponding CouenneObject was created on w=f(x,y), it is
  /// either x or y, chosen previously with a call to getFixVar()
  /// expression *reference_;
  expression *variable_;

  /// SmartPointer to the Journalist
  JnlstPtr jnlst_;

  /// shall we do Feasibility based Bound Tightening (FBBT) at branching?
  bool doFBBT_;

  /// shall we add convexification cuts at branching?
  bool doConvCuts_;

  /// down branch estimate (done at selectBranch with reduced costs)
  double downEstimate_;

  /// up   branch estimate
  double upEstimate_;
};

#endif

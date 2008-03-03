/*
 * Name:    CouenneObject.hpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Object for auxiliary variables
 *
 * (C) Carnegie-Mellon University, 2006-07.
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNEOBJECT_HPP
#define COUENNEOBJECT_HPP

#include "BonBabSetupBase.hpp"
#include "CoinFinite.hpp"
#include "OsiBranchingObject.hpp"

#include "exprVar.hpp"

#include "CouenneJournalist.hpp"

#define AGGR_MUL 2
const CouNumber closeToBounds = .01;


/// Define what kind of branching (two- or three-way) and where to
/// start from: left, (center,) or right. The last is to help
/// diversify branching through randomization, which may help when the
/// same variable is branched upon in several points of the BB tree.

enum {TWO_LEFT,                 TWO_RIGHT,   TWO_RAND,
      THREE_LEFT, THREE_CENTER, THREE_RIGHT, THREE_RAND, BRANCH_NONE};

//
class funtriplet;
CouNumber minMaxDelta (funtriplet *ft, CouNumber lb, CouNumber ub);
CouNumber maxHeight   (funtriplet *ft, CouNumber lb, CouNumber ub);

/// OsiObject for auxiliary variables $w=f(x)$. 
///
/// Associated with a multi-variate function $f(x)$ and a related
/// infeasibility $|w-f(x)|$, creates branches to help restoring
/// feasibility

class CouenneObject: public OsiObject {

public:

  /// strategy names
  enum brSelStrat {NO_BRANCH, MID_INTERVAL, MIN_AREA, BALANCED, LP_CENTRAL, LP_CLAMPED};

  /// Constructor with information for branching point selection strategy
  CouenneObject (exprVar *ref, Bonmin::BabSetupBase *base,
		 JnlstPtr jnlst);

  /// Destructor
  ~CouenneObject () {}

  /// Copy constructor
  CouenneObject (const CouenneObject &src);

  /// Cloning method
  virtual OsiObject * clone () const
  {return new CouenneObject (*this);}

  /// compute infeasibility of this variable, |w - f(x)| (where w is
  /// the auxiliary variable defined as w = f(x)
  /// TODO: suggest way
  virtual inline double infeasibility (const OsiBranchingInformation *info, int &way) const {  

    if (strategy_ == NO_BRANCH) return 0.;

    CouNumber delta = 
      fabs (info -> solution_ [reference_ -> Index ()] - 
	    (*(reference_ -> Image ())) ());

    return (delta < CoinMin (COUENNE_EPS, feas_tolerance_)) ? 0. : delta;
  }

  /// fix (one of the) arguments of reference auxiliary variable 
  virtual double feasibleRegion (OsiSolverInterface*, 
				 const OsiBranchingInformation*) const;

  /// create CouenneBranchingObject or CouenneThreeWayBranchObj based
  /// on this object
  virtual OsiBranchingObject* createBranch (OsiSolverInterface*, 
					    const OsiBranchingInformation*, int) const;

  /// return reference auxiliary variable
  exprVar *Reference () const
  {return reference_;}

  /// return branching point selection strategy
  enum brSelStrat Strategy ()  const
  {return strategy_;}

  /// pick branching point based on current strategy
  CouNumber getBrPoint (funtriplet *ft, CouNumber x0, CouNumber l, CouNumber u) const;

  /// returns a point "inside enough" a given interval, or x if it is
  /// already
  CouNumber midInterval (CouNumber x, CouNumber l, CouNumber u) const;

protected:

  /// The (auxiliary) variable this branching object refers to. If the
  /// expression is w=f(x,y), this is w, as opposed to
  /// CouenneBranchingObject, where it would be either x or y.
  exprVar *reference_;

  /// Branching point selection strategy
  enum brSelStrat strategy_;

  /// SmartPointer to the Journalist
  JnlstPtr jnlst_;

  /// Combination parameter for the mid-point branching point
  /// selection strategy
  CouNumber alpha_;

  /// Defines safe interval percentage for using LP point as a branching point
  CouNumber lp_clamp_;

  /// feasibility tolerance (equal to that of CouenneProblem)
  CouNumber feas_tolerance_;

  /// shall we do Feasibility based Bound Tightening (FBBT) at branching?
  bool doFBBT_;

  /// shall we add convexification cuts at branching?
  bool doConvCuts_;
};

#endif

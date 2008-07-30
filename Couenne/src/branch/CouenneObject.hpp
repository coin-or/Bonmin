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

#include "exprVar.hpp"
#include "CouenneJournalist.hpp"
#include "OsiBranchingObject.hpp"

#define AGGR_MUL 2
#define THRES_ZERO_SYMM 0.8
#define TOL 0.

const CouNumber closeToBounds = .05;

/// Define what kind of branching (two- or three-way) and where to
/// start from: left, (center,) or right. The last is to help
/// diversify branching through randomization, which may help when the
/// same variable is branched upon in several points of the BB tree.

enum {TWO_LEFT,                 TWO_RIGHT,   TWO_RAND,
      THREE_LEFT, THREE_CENTER, THREE_RIGHT, THREE_RAND, BRANCH_NONE};

class funtriplet;
class CouenneProblem;

CouNumber minMaxDelta (funtriplet *ft, CouNumber lb, CouNumber ub);
CouNumber maxHeight   (funtriplet *ft, CouNumber lb, CouNumber ub);


/// OsiObject for auxiliary variables $w=f(x)$. 
///
/// Associated with a multi-variate function $f(x)$ and a related
/// infeasibility $|w-f(x)|$, creates branches to help restoring
/// feasibility

class CouenneObject: public OsiObject {

public:

  /// type of up/down estimate to return for pseudocosts
  enum pseudocostMult {INFEASIBILITY, 
		       INTERVAL_LP, INTERVAL_LP_REV,
		       INTERVAL_BR, INTERVAL_BR_REV,
		       PROJECTDIST};

  /// type of object (for branching variable selection)
  enum branch_obj {EXPR_OBJ, VAR_OBJ, VT_OBJ};

  /// strategy names
  enum brSelStrat {NO_STRATEGY, NO_BRANCH, MID_INTERVAL, MIN_AREA, BALANCED, LP_CENTRAL, LP_CLAMPED};

  /// empty constructor (for unused objects)
  CouenneObject ();

  /// Constructor with information for branching point selection strategy
  CouenneObject (CouenneProblem *p, 
		 exprVar *ref, Bonmin::BabSetupBase *base, JnlstPtr jnlst);

  /// Constructor with lesser information, used for infeasibility only
  CouenneObject (exprVar *ref, Bonmin::BabSetupBase *base, JnlstPtr jnlst);

  /// Destructor
  ~CouenneObject () {}

  /// Copy constructor
  CouenneObject (const CouenneObject &src);

  /// Cloning method
  virtual OsiObject * clone () const
  {return new CouenneObject (*this);}

  /// set object parameters by reading from command line
  void setParameters (Bonmin::BabSetupBase *base);

  /// compute infeasibility of this variable, |w - f(x)| (where w is
  /// the auxiliary variable defined as w = f(x)
  virtual double infeasibility (const OsiBranchingInformation *info, int &way) const;

  /// compute infeasibility of this variable, |w - f(x)|, where w is
  /// the auxiliary variable defined as w = f(x)
  virtual double checkInfeasibility (const OsiBranchingInformation * info) const;

  /// fix (one of the) arguments of reference auxiliary variable 
  virtual double feasibleRegion (OsiSolverInterface*, const OsiBranchingInformation*) const;

  /// create CouenneBranchingObject or CouenneThreeWayBranchObj based
  /// on this object
  virtual OsiBranchingObject *createBranch (OsiSolverInterface*, 
					    const OsiBranchingInformation*, int) const;

  /// return reference auxiliary variable
  exprVar *Reference () const
  {return reference_;}

  /// return branching point selection strategy
  enum brSelStrat Strategy ()  const
  {return strategy_;}

  /// pick branching point based on current strategy
  CouNumber getBrPoint (funtriplet *ft, CouNumber x0, CouNumber l, CouNumber u) const;

  /// returns a point "inside enough" a given interval, or x if it
  /// already is
  CouNumber midInterval (CouNumber x, CouNumber l, CouNumber u) const;

  /// Return "down" estimate (for non-convex, distance old <--> new LP point)
  virtual double downEstimate () const
  {//if (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING)) {
    //printf ("DOWN EST = %g for ", downEstimate_); 
    //reference_ -> print (); 
    //printf ("\n");
    //}
  return downEstimate_;}

  /// Return "up" estimate (for non-convex, distance old <--> new LP point)
  virtual double upEstimate () const
  {//if (jnlst_ -> ProduceOutput (J_MATRIX, J_BRANCHING)) {
    //printf ("UP EST = %g for ", upEstimate_); 
    //reference_ -> print (); 
    //printf ("\n");
    //}
  return upEstimate_;}

  /// set up/down estimate (0 for down, 1 for up). This happens in
  /// CouenneChooseStrong, where a new LP point is available and we
  /// can measure distance from old LP point. This is the denominator
  /// we use in pseudocost
  void setEstimate (double est, int direction)
  {(direction ? upEstimate_ : downEstimate_) = est;}

  /// set up/down estimates based on branching information
  void setEstimates (const OsiBranchingInformation *info,
		     CouNumber *infeasibility,
		     CouNumber *brpt) const;

  /// are we on the bad or good side of the expression?
  virtual bool isCuttable () const
  {return reference_ -> Image () -> isCuttable (problem_, reference_ -> Index ());}

protected:

  /// pointer to Couenne problem
  CouenneProblem *problem_;

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

  /// down estimate (to be used in pseudocost)
  mutable double downEstimate_;

  /// up estimate (to be used in pseudocost)
  mutable double upEstimate_;

  /// multiplier type for pseudocost
  enum pseudocostMult pseudoMultType_;
};

#endif

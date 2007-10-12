/*
 * Name:    CouenneObject.hpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Object for auxiliary variables
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNEOBJECT_HPP
#define COUENNEOBJECT_HPP

#include <CoinFinite.hpp>
#include <OsiBranchingObject.hpp>
#include <exprAux.hpp>


/// Define what kind of branching (two- or three-way) and where to
/// start from: left, (center,) or right. The last is to help diversify
/// branching through randomization, which may help when the same
/// variable is branched upon in several points of the BB tree.

enum {TWO_LEFT,                 TWO_RIGHT,   TWO_RAND,
      THREE_LEFT, THREE_CENTER, THREE_RIGHT, THREE_RAND, BRANCH_NONE};


/// OsiObject for auxiliary variables $w=f(x)$. 
///
/// Associated with a multi-variate function $f(x)$ and a related
/// infeasibility $|w=f(x)|$, creates branches to help restoring
/// feasibility

class CouenneObject: public OsiObject {

public:

  /// Constructor
  CouenneObject (exprAux *ref):
    reference_ (ref),
    brVarInd_  (-1), 
    brPts_     (NULL),
    whichWay_  (BRANCH_NONE) {}

  /// Constructor
  CouenneObject (exprVar *ref):
    reference_ (NULL),
    brVarInd_  (-1), 
    brPts_     (NULL),
    whichWay_  (BRANCH_NONE) 
  {printf ("invalid CouenneObject: original variable %d\n", ref -> Index ()); exit (-1);}

  /// Destructor
  ~CouenneObject () 
  {if (brPts_) free (brPts_);}

  /// Copy constructor
  CouenneObject (CouenneObject &src):
    reference_ (src.reference_),
    brVarInd_  (src.brVarInd_),
    brPts_     (NULL) {

    if (src.brPts_) {
      register int i=0;
      while (src.brPts_ [i] > -COIN_DBL_MAX) i++;

      if (i) {
	brPts_ = (CouNumber *) malloc (i * sizeof (CouNumber));
	while (i--) 
	  brPts_ [i] = src.brPts_ [i];
      }
    }
  }

  /// Cloning method
  virtual OsiObject * clone() const
  {return new CouenneObject (reference_);}

  /// compute infeasibility of this variable, |w - f(x)| (where w is
  /// the auxiliary variable defined as w = f(x)
  virtual double infeasibility (const OsiBranchingInformation*, int &) const;

  /// fix (one of the) arguments of reference auxiliary variable 
  virtual double feasibleRegion (OsiSolverInterface*, 
				 const OsiBranchingInformation*) const;

  /// create CouenneBranchingObject or CouenneThreeWayBranchObj based
  /// on this object
  virtual OsiBranchingObject* createBranch (OsiSolverInterface*, 
					    const OsiBranchingInformation*, int) const;

  /// return reference auxiliary variable
  exprAux *Reference ()
  {return reference_;}

protected:

  /// The (auxiliary) variable which this branching object refers
  /// to. If the expression is w=f(x,y), this is w, as opposed to
  /// CouenneBranchingObject, where it would be either x or y.
  exprAux *reference_;

  /// index on the branching variable
  mutable int brVarInd_;

  /// where to branch. It is a vector in the event we want to use a
  /// ThreeWayBranching. Ends with a -COIN_DBL_MAX (not considered...)
  mutable CouNumber *brPts_;

  /// How to branch. This should be done automatically from
  /// ::infeasibility() and ::createBranch(), but for some reason
  /// obj->whichWay() in the last 20 lines of CbcNode.cpp returns 0,
  /// therefore we set our own whichWay_ within this object and use it
  /// between the two methods.
  mutable int whichWay_;
};

#endif

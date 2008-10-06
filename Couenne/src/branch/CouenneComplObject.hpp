/*
 * Name:    CouenneComplObject.hpp
 * Authors: Pietro Belotti, Lehigh University
 * Purpose: Branching object for complementarity constraints
 *
 * (C) Carnegie-Mellon University, 2008.
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNECOMPLOBJECT_HPP
#define COUENNECOMPLOBJECT_HPP

#include "CouenneObject.hpp"

/// OsiObject for complementarity constraints $x_1 x_2 = 0$. 
///
/// Associated with two variables x_1 and x_2, branches with either x_1=0 or x_2=0

class CouenneComplObject: public CouenneObject {

public:

  /// empty constructor (for unused objects)
  CouenneComplObject ();

  /// Constructor with information for branching point selection strategy
  CouenneComplObject (CouenneProblem *p, 
		      exprVar *ref, Bonmin::BabSetupBase *base, JnlstPtr jnlst);

  /// Constructor with lesser information, used for infeasibility only
  CouenneComplObject (exprVar *ref, Bonmin::BabSetupBase *base, JnlstPtr jnlst);

  /// Destructor
  ~CouenneComplObject () {}

  /// Copy constructor
  CouenneComplObject (const CouenneObject &src);

  /// Cloning method
  virtual OsiObject * clone () const
  {return new CouenneComplObject (*this);}

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

  /// are we on the bad or good side of the expression?
  virtual bool isCuttable () const {
    return (reference_ -> Image ()) ? 
      ((!(reference_ -> isInteger ())) &&
       reference_ -> Image () -> isCuttable (problem_, reference_ -> Index ())) :
      (!(reference_ -> isInteger ()));
  }

protected:

  /// The (auxiliary) variable this branching object refers to. If the
  /// expression is w=f(x,y), this is w, as opposed to
  /// CouenneBranchingObject, where it would be either x or y.
  exprVar *reference_;
};

#endif

/*
 * Name:    CouenneVarObject.hpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: Object for branching on variables
 *
 * (C) Carnegie-Mellon University, 2008.
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNEVAROBJECT_HPP
#define COUENNEVAROBJECT_HPP

#include "CouenneObject.hpp"
#include "CouenneProblem.hpp"


/// OsiObject for variables in a MINLP
class CouenneVarObject: public CouenneObject {

public:

  /// Constructor with information for branching point selection strategy
  CouenneVarObject (exprVar *ref, 
		    CouenneProblem *p,
		    Bonmin::BabSetupBase *base, 
		    JnlstPtr jnlst);

  /// Copy constructor
  CouenneVarObject (const CouenneVarObject &src);

  /// Destructor
  ~CouenneVarObject () {}

  /// Cloning method
  virtual OsiObject * clone () const
  {return new CouenneVarObject (*this);}

  /// compute infeasibility of this variable x as the sum/min/max of
  /// all infeasibilities of auxiliaries w whose defining function
  /// depends on x |w - f(x)|
  ///
  /// TODO: suggest way
  virtual double infeasibility (const OsiBranchingInformation *info, int &way) const;

  /// create CouenneBranchingObject or CouenneThreeWayBranchObj based
  /// on this object
  virtual OsiBranchingObject* createBranch (OsiSolverInterface*, 
					    const OsiBranchingInformation*, int) const;

  /// Return "down" estimate (for non-convex, distance old <--> new LP point)
  virtual double downEstimate () const
  {//printf("downest = %e\n", downEstimate_); 
    return downEstimate_;}

  /// Return "up" estimate (for non-convex, distance old <--> new LP point)
  virtual double upEstimate () const 
  {//printf("upest = %e\n", upEstimate_); 
    return upEstimate_;}

  /// set up/down estimate (0 for down, 1 for up). This happens in
  /// CouenneChooseStrong, where a new LP point is available and we
  /// can measure distance from old LP point. This is the denominator
  /// we use in pseudocost
  void setEstimate (double est, int direction)
  {(direction ? upEstimate_ : downEstimate_) = est;}

  void setProblem (CouenneProblem *p)
  {problem_ = p;}

protected:

  /// up estimate (to be used in pseudocost)
  mutable double upEstimate_;

  /// down estimate (to be used in pseudocost)
  mutable double downEstimate_;

  /// pointer to problem
  CouenneProblem *problem_;

  /// Method computing the branching point
  CouNumber computeBranchingPoint(const OsiBranchingInformation *info,
				 int& bestWay) const;
};

#endif

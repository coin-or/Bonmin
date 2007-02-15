/*
 * Name:    CouenneObject.hpp
 * Author:  Pietro Belotti
 * Purpose: Object for auxiliary variables
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNEOBJECT_HPP
#define COUENNEOBJECT_HPP

#include <CoinFinite.hpp>
#include <OsiBranchingObject.hpp>
#include <exprAux.h>


class CouenneObject: public OsiObject {

public:

  /// Constructor
  CouenneObject (exprAux *ref):
    reference_ (ref) {}

  /// Cloning method
  virtual OsiObject * clone() const
  {return new CouenneObject (reference_);}

  /// compute infeasibility of this variable, |w - f(x)| (where w is
  /// the auxiliary variable defined as w = f(x)
  virtual double infeasibility (const OsiBranchingInformation*, int &) const;

  /// fix (one of the) argument of reference auxiliary variable 
  virtual double feasibleRegion (OsiSolverInterface*, const OsiBranchingInformation*) const;

  /// create CouenneBranchingObject based on this object
  virtual OsiBranchingObject* createBranch (OsiSolverInterface*, 
					    const OsiBranchingInformation*, int) const;

protected:

  /// the (auxiliary) variable which this branching object refers to
  exprAux *reference_;
};

#endif

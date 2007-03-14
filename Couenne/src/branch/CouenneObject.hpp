/*
 * Name:    CouenneObject.hpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
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

  /// Destructor
  ~CouenneObject () {}

  /// Copy constructor
  CouenneObject (CouenneObject &src):
    reference_ (src.reference_) {}

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

  /// return reference auxiliary variable
  exprAux *Reference () 
  {return reference_;}

protected:

  /// The (auxiliary) variable which this branching object refers
  /// to. If the expression is w=f(x,y), this is w, as opposed to
  /// CouenneBranchingObject, where it would be either x or y.
  exprAux *reference_;
};

#endif

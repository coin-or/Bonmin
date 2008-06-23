/*
 * Name:    CouenneVTObject.hpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: Object for branching on variables using violation transfer
 *
 * (C) Carnegie-Mellon University, 2008.
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNEVTOBJECT_HPP
#define COUENNEVTOBJECT_HPP

#include "CouenneVarObject.hpp"

/// OsiObject for violation transfer on variables in a MINLP
class CouenneVTObject: public CouenneVarObject {

public:

  /// Constructor with information for branching point selection strategy
  CouenneVTObject (CouenneProblem *p,
		   exprVar *ref, 
		   Bonmin::BabSetupBase *base, 
		   JnlstPtr jnlst):

    CouenneVarObject (p, ref, base, jnlst) {}

  /// Copy constructor
  CouenneVTObject (const CouenneVTObject &src):
    CouenneVarObject (src) {}

  /// Destructor
  ~CouenneVTObject () {}

  /// Cloning method
  virtual OsiObject *clone () const
  {return new CouenneVTObject (*this);}

  /// compute infeasibility of this variable x as the sum/min/max of
  /// all infeasibilities of auxiliaries w whose defining function
  /// depends on x |w - f(x)|
  virtual double infeasibility (const OsiBranchingInformation *info, int &way) const;
};

#endif

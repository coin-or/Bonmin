// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
//
// Author:  Andreas Waechter      2007-08-20    IBM
//
#ifndef BonStrongBranchingSolver_H
#define BonStrongBranchingSolver_H

#include "BonOsiTMINLPInterface.hpp"
#include "BonRegisteredOptions.hpp"
namespace Bonmin {

/** This class is the base class for a solver that can be used in
 *  BonOsiSolverInterface to perform the strong branching solves.
*/

class StrongBranchingSolver : public ReferencedObject  {
 
public:

  /// Constructor from solver
  StrongBranchingSolver (OsiTMINLPInterface * solver);

  /// Assignment operator 
  StrongBranchingSolver & operator= (const StrongBranchingSolver& rhs);
  /// Copy constructor
  StrongBranchingSolver(const StrongBranchingSolver& rhs);

  /// Destructor
  virtual ~StrongBranchingSolver ();

  /// Called to initialize solver before a bunch of strong branching
  /// solves
  virtual void markHotStart(OsiTMINLPInterface* tminlp_interface) = 0;

  /// Called to solve the current TMINLP (with changed bound information)
  virtual TNLPSolver::ReturnStatus solveFromHotStart(OsiTMINLPInterface* tminlp_interface) = 0;

  /// Called after all strong branching solves in a node
  virtual void unmarkHotStart(OsiTMINLPInterface* tminlp_interface) = 0;

protected:

  inline SmartPtr<Journalist>& Jnlst()
  {
    return jnlst_;
  }
  inline SmartPtr<OptionsList>& Options()
  {
    return options_;
  }
  inline SmartPtr<RegisteredOptions>& RegOptions()
  {
    return reg_options_;
  }
private:
  /** Default Constructor, forbiden for some reason.*/
  StrongBranchingSolver ();

  SmartPtr<Journalist> jnlst_;
  SmartPtr<OptionsList> options_;
  SmartPtr<Bonmin::RegisteredOptions> reg_options_;

  int bb_log_level_;

};

}
#endif

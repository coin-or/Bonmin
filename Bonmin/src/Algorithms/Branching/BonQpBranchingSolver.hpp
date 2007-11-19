// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef BonQpBranchingSolver_H
#define BonQpBranchingSolver_H

#include "BonStrongBranchingSolver.hpp"
#include "BonBranchingTQP.hpp"

namespace Bonmin
{

  /** This class chooses a variable to branch on

      This implementation solves the Qp model for different branches
      (strong branching).
  */

  class QpBranchingSolver : public StrongBranchingSolver
  {

  public:

    /// Constructor from solver (so we can set up arrays etc)
    QpBranchingSolver (OsiTMINLPInterface * solver);

    /// Copy constructor
    QpBranchingSolver (const QpBranchingSolver &);

    /// Assignment operator
    QpBranchingSolver & operator= (const QpBranchingSolver& rhs);

    /// Destructor
    virtual ~QpBranchingSolver ();

    /// Called to initialize solver before a bunch of strong branching
    /// solves
    virtual void markHotStart(OsiTMINLPInterface* tminlp_interface);

    /// Called to solve the current TMINLP (with changed bound information)
    virtual TNLPSolver::ReturnStatus solveFromHotStart(OsiTMINLPInterface* tminlp_interface);

    /// Called after all strong branching solves in a node
    virtual void unmarkHotStart(OsiTMINLPInterface* tminlp_interface);

  private:
    /// Default Constructor
    QpBranchingSolver ();

    SmartPtr<BranchingTQP> branching_tqp_;

    SmartPtr<TNLPSolver> tqp_solver_;

    bool first_solve_;
  };

}

#endif

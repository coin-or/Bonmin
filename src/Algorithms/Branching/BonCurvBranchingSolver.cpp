// Copyright (C) 2006, 2007 International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "BonCurvBranchingSolver.hpp"

namespace Bonmin
{

  CurvBranchingSolver::CurvBranchingSolver(OsiTMINLPInterface * solver) :
      StrongBranchingSolver(solver),
      orig_d_(NULL),
      projected_d_(NULL),
      x_l_orig_(NULL),
      x_u_orig_(NULL),
      g_l_orig_(NULL),
      g_u_orig_(NULL),
      solution_(NULL),
      duals_(NULL)
  {}

  CurvBranchingSolver::CurvBranchingSolver(const CurvBranchingSolver & rhs) :
      StrongBranchingSolver(rhs),
      orig_d_(NULL),
      projected_d_(NULL),
      x_l_orig_(NULL),
      x_u_orig_(NULL),
      g_l_orig_(NULL),
      g_u_orig_(NULL),
      solution_(NULL),
      duals_(NULL)
  {}

  CurvBranchingSolver &
  CurvBranchingSolver::operator=(const CurvBranchingSolver & rhs)
  {
    assert(!x_l_orig_);
    if (this != &rhs) {
      StrongBranchingSolver::operator=(rhs);
    }
    return *this;
  }

  CurvBranchingSolver::~CurvBranchingSolver ()
  {
    delete [] orig_d_;
    delete [] projected_d_;
    delete [] x_l_orig_;
    delete [] x_u_orig_;
    delete [] g_l_orig_;
    delete [] g_u_orig_;
    delete [] solution_;
    delete [] duals_;
  }

  void CurvBranchingSolver::
  markHotStart(OsiTMINLPInterface* tminlp_interface)
  {
    if (IsNull(cur_estimator_)) {
      // Get a curvature estimator
      cur_estimator_ = new CurvatureEstimator(Jnlst(), Options(),
          tminlp_interface->problem());
    }

    new_bounds_ = true;
    new_x_ = true;
    new_mults_ = true;

    delete [] solution_;
    delete [] duals_;
    solution_ =  NULL;
    duals_ =  NULL;

    numCols_ = tminlp_interface->getNumCols();
    numRows_ = tminlp_interface->getNumRows();
    solution_ = CoinCopyOfArray(tminlp_interface->problem()->x_sol(), numCols_);
    duals_ = CoinCopyOfArray(tminlp_interface->problem()->duals_sol(),
        numRows_ + 2*numCols_);
    obj_value_ = tminlp_interface->problem()->obj_value();

    delete [] orig_d_;
    delete [] projected_d_;
    orig_d_ = NULL;
    projected_d_ = NULL;
    orig_d_ = new double[numCols_];
    projected_d_ = new double[numCols_];

    // Get a copy of the current bounds
    delete [] x_l_orig_;
    delete [] x_u_orig_;
    delete [] g_l_orig_;
    delete [] g_u_orig_;
    x_l_orig_ = NULL;
    x_u_orig_ = NULL;
    g_l_orig_ = NULL;
    g_u_orig_ = NULL;
    x_l_orig_ = new Number[numCols_];
    x_u_orig_ = new Number[numCols_];
    g_l_orig_ = new Number[numRows_];
    g_u_orig_ = new Number[numRows_];

#ifndef NDEBUG
    bool retval =
#endif
      tminlp_interface->problem()->
      get_bounds_info(numCols_, x_l_orig_, x_u_orig_,
          numRows_, g_l_orig_, g_u_orig_);
    assert(retval);
  }

  void CurvBranchingSolver::
  unmarkHotStart(OsiTMINLPInterface* tminlp_interface)
  {
    // Free memory
    delete [] solution_;
    delete [] duals_;
    solution_ =  NULL;
    duals_ =  NULL;
    delete [] orig_d_;
    delete [] projected_d_;
    orig_d_ = NULL;
    projected_d_ = NULL;
    delete [] x_l_orig_;
    delete [] x_u_orig_;
    delete [] g_l_orig_;
    delete [] g_u_orig_;
    x_l_orig_ = NULL;
    x_u_orig_ = NULL;
    g_l_orig_ = NULL;
    g_u_orig_ = NULL;
  }

  TNLPSolver::ReturnStatus CurvBranchingSolver::
  solveFromHotStart(OsiTMINLPInterface* tminlp_interface)
  {
    // return iteration limit reached as status, so that it is
    // clear we don't have a feasible point
    TNLPSolver::ReturnStatus retstatus = TNLPSolver::iterationLimit;

    const double* z_L = duals_;
    const double* z_U = z_L + numCols_;
    const double* lam = z_U + numCols_;

    // get the current points
    const double* b_L = tminlp_interface->getColLower();
    const double* b_U = tminlp_interface->getColUpper();

    // Compute the step from the current point to the solution
    // ToDo: If we know what changes, we can be more efficient
    for (int i=0; i<numCols_; i++) {
      if (b_L[i]>solution_[i]) {
        orig_d_[i] = solution_[i]-b_L[i];
      }
      else if (b_U[i]<solution_[i]) {
        orig_d_[i] = solution_[i]-b_U[i];
      }
      else {
        orig_d_[i] = 0.;
      }
    }

    double gradLagTd;
    double dTHLagd;
    bool retval =
      cur_estimator_->ComputeNullSpaceCurvature(
        numCols_, solution_, new_x_, x_l_orig_, x_u_orig_,
        g_l_orig_, g_u_orig_, new_bounds_, z_L, z_U,
        numRows_, lam, new_mults_, orig_d_, projected_d_,
        gradLagTd, dTHLagd);

    if (!retval) {
      retstatus = TNLPSolver::computationError;
    }
    else {
      new_bounds_ = false;
      new_x_ = false;
      new_mults_ = false;
      const double alpha = 1.0; // Think about this
      double new_obj_value = obj_value_ +
          alpha*gradLagTd + 0.5*alpha*alpha*dTHLagd;
      tminlp_interface->problem()->set_obj_value(new_obj_value);
    }

    return retstatus;
  }


}

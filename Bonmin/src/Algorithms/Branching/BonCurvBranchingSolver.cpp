// Copyright (C) 2006, 2007 International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "BonCurvBranchingSolver.hpp"

namespace Bonmin {

CurvBranchingSolver::CurvBranchingSolver(OsiTMINLPInterface * solver) :
  StrongBranchingSolver(solver),
  orig_d_(NULL),
  projected_d_(NULL)
{
}

CurvBranchingSolver::CurvBranchingSolver(const CurvBranchingSolver & rhs) :
  StrongBranchingSolver(rhs),
  orig_d_(NULL),
  projected_d_(NULL)
{
}

CurvBranchingSolver &
CurvBranchingSolver::operator=(const CurvBranchingSolver & rhs)
{
  if (this != &rhs) {
    StrongBranchingSolver::operator=(rhs);
  }
  return *this;
}

CurvBranchingSolver::~CurvBranchingSolver ()
{
  delete [] orig_d_;
  delete [] projected_d_;
}

void CurvBranchingSolver::
markHotStart(OsiTMINLPInterface* tminlp_interface)
{
  // Get a new curvature estimator (or should we just keep one?)
  cur_estimator_ = new CurvatureEstimator(Jnlst(), Options(),
					  tminlp_interface->problem());
  new_bounds_ = true;
  new_x_ = true;
  new_mults_ = true;

  numCols_ = tminlp_interface->getNumCols();
  numRows_ = tminlp_interface->getNumRows();
  solution_ = tminlp_interface->problem()->x_sol();
  duals_ = tminlp_interface->problem()->duals_sol();
  obj_value_ = tminlp_interface->problem()->obj_value();

  delete [] orig_d_;
  delete [] projected_d_;
  orig_d_ = NULL;
  projected_d_ = NULL;
  orig_d_ = new double[numCols_];
  projected_d_ = new double[numCols_];
}

void CurvBranchingSolver::
unmarkHotStart(OsiTMINLPInterface* tminlp_interface)
{
  // Free memory
  cur_estimator_ = NULL;
  delete [] orig_d_;
  delete [] projected_d_;
  orig_d_ = NULL;
  projected_d_ = NULL;
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
    orig_d_[i] = max(max(0., b_L[i]-solution_[i]), solution_[i]-b_U[i]);
  }

  double gradLagTd;
  double dTHLagd;
  bool retval =
    cur_estimator_->ComputeNullSpaceCurvature(
           new_bounds_, numCols_, solution_, new_x_, z_L, z_U,
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

// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "BonCurvBranching.hpp"
#include "IpBlas.hpp"

namespace Bonmin {

BonCurvBranching::BonCurvBranching(OsiTMINLPInterface * solver) :
  BonChooseVariable(solver)
{
  SmartPtr<TNLPSolver> tnlp_solver =
    dynamic_cast<TNLPSolver *> (solver->solver());
  SmartPtr<OptionsList> options = tnlp_solver->Options();
  SmartPtr<TNLP> tnlp = solver->problem();

  cur_estimator_ = new CurvatureEstimator(jnlst_, options, tnlp);
}

BonCurvBranching::BonCurvBranching(const BonCurvBranching & rhs) :
  BonChooseVariable(rhs)
{
  cur_estimator_ = rhs.cur_estimator_;
}

BonCurvBranching &
BonCurvBranching::operator=(const BonCurvBranching & rhs)
{
  if (this != &rhs) {
    BonChooseVariable::operator=(rhs);
    cur_estimator_ = rhs.cur_estimator_;
  }
  return *this;
}

BonCurvBranching::~BonCurvBranching ()
{}

// Clone
OsiChooseVariable *
BonCurvBranching::clone() const
{
  return new BonCurvBranching(*this);
}

int
BonCurvBranching::fill_changes(OsiSolverInterface * solver,
			       OsiBranchingInformation *info,
			       bool fixVariables, int numStrong,
			       double* change_down,
			       double* change_up, int& best_way)
{
  // Get info about the current solution
  const double* solution = solver->getColSolution();// Current solution
  int numCols = solver->getNumCols();
  int numRows = solver->getNumRows();
  const double* lam = solver->getRowPrice();
  const double* z_L = lam + numRows;
  const double* z_U = z_L + numCols;
  const double* b_L = solver->getColLower();
  const double* b_U = solver->getColUpper();

  // Set up stuff for the curvature estimator
  bool new_bounds = true;
  bool new_x = true;
  bool new_mults = true;
  double * orig_d = new double[numCols];
  double * projected_d = new double[numCols];
  const Number zero = 0.;
  IpBlasDcopy(numCols, &zero, 0, orig_d, 1);

  const Number large_number = COIN_DBL_MAX;
  for (int i=0; i<numStrong; i++) {
    int index = list_[i];
    const OsiObject * object = solver->object(index);
    int col_number = object->columnNumber();
    DBG_ASSERT(col_number != -1);
    // up
    orig_d[col_number] = 1.;
    double gradLagTd;
    double dTHLagd;
    bool retval =
      cur_estimator_->ComputeNullSpaceCurvature(
         new_bounds, numCols, solution, new_x, z_L, z_U,
	 numRows, lam, new_mults, orig_d, projected_d, gradLagTd, dTHLagd);
    // ToDo
    if (!retval) {
      printf("Problem with curvature estimator for down up.\n");
    }
    new_bounds = false;
    new_x = false;
    new_mults = false;

    // Determine step size and predicted change
    const double &curr_val = solution[col_number];
    if (retval && projected_d[col_number] != 0.) {
      const double up_val = Min(b_U[col_number],ceil(curr_val));
      double alpha = (up_val-curr_val)/projected_d[col_number];
      change_up[i] = alpha*gradLagTd + 0.5*alpha*alpha*dTHLagd;
    }
    else {
      change_up[i] = -large_number;
    }

    // down
    orig_d[col_number] = -1.;
    retval = cur_estimator_->ComputeNullSpaceCurvature(
          new_bounds, numCols, solution, new_x, z_L, z_U,
	  numRows, lam, new_mults, orig_d, projected_d, gradLagTd, dTHLagd);
    // ToDo
    if (!retval) {
      printf("Problem with curvature estimator for down down.\n");
    }

    // Determine step size and predicted change
    if (retval && projected_d[col_number] != 0.) {
      const double down_val = Max(b_L[col_number],floor(curr_val));
      double alpha = (down_val-curr_val)/projected_d[col_number];
      change_down[i] = alpha*gradLagTd + 0.5*alpha*alpha*dTHLagd;
    }
    else {
      change_down[i] = -large_number;
    }

    orig_d[col_number] = 0.;
  }

  delete [] orig_d;
  delete [] projected_d;

  return -1; // We never detect an infeasible problem
}

}

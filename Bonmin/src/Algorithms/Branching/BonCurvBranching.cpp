// Copyright (C) 2006, 2007 International Business Machines
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

// Clone
ChooseVariable *
BonCurvBranching::clone2() const
{
  return new BonCurvBranching(*this);
}

static int find_where_to_branch_sos(
  const double* solution,
  OsiBranchingInformation *info,
  const OsiSOS* sos_object
  )
{
  // This code is copied from OsiSOS::createBranch
  int j;
  double tolerance = 1e-8; //TODO: What goes here?
  const double * upper = info->upper_;
  int firstNonFixed=-1;
  int lastNonFixed=-1;
  int firstNonZero=-1;
  int lastNonZero = -1;
  double weight = 0.0;
  double sum =0.0;
  const int* members = sos_object->members();
  const double* weights = sos_object->weights();
  for (j=0;j<sos_object->numberMembers();j++) {
    int iColumn = members[j];
    if (upper[iColumn]) {
      double value = CoinMax(0.0,solution[iColumn]);
      sum += value;
      if (firstNonFixed<0)
	firstNonFixed=j;
      lastNonFixed=j;
      if (value>tolerance) {
	weight += weights[j]*value;
	if (firstNonZero<0)
	  firstNonZero=j;
	lastNonZero=j;
      }
    }
  }
  // find where to branch
  assert (sum>0.0);
  weight /= sum;
  int iWhere;
  for (iWhere=firstNonZero;iWhere<lastNonZero;iWhere++) 
    if (weight<weights[iWhere+1])
      break;
  return iWhere;
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
    if (col_number != -1) {
      // This is a regular integer variable
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
	printf("Problem with curvature estimator for up.\n");
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
	printf("Problem with curvature estimator for down.\n");
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
    else {
      // SOS Branching object
      const OsiSOS* sos_object
	= dynamic_cast<const OsiSOS*>(object);
      if (sos_object->sosType() != 1) {
	printf("Currently, curvature branching can only handle type-1 SOS constraints\n");
	exit(-1);
      }
      const int* members = sos_object->members();
      //for(int j=0; j<sos_object->numberMembers(); j++) {
      //printf("%3d %25.16e\n", members[j], solution[members[j]]);
      //}
      int iWhere = find_where_to_branch_sos(solution, info, sos_object);
      //printf("iWhere = %d \n", iWhere);
      // up (fix all members[i] with i<=iWhere)
      const int nMembers = sos_object->numberMembers();
      for (int j=0; j<=iWhere; j++) {
	const int& iCol = members[j];
	orig_d[iCol] = -solution[iCol];
      }
      double sum = 0.;
      for (int j=iWhere+1; j<nMembers; j++) {
	sum += solution[members[j]];
      }
      double d = (1.-sum)/((double)(nMembers-iWhere-1));
      for (int j=iWhere+1; j<nMembers; j++) {
	orig_d[members[j]] = d;
      }

      double gradLagTd;
      double dTHLagd;
      bool retval =
	cur_estimator_->ComputeNullSpaceCurvature(
           new_bounds, numCols, solution, new_x, z_L, z_U,
	   numRows, lam, new_mults, orig_d, projected_d, gradLagTd, dTHLagd);
      // ToDo
      if (!retval) {
	printf("Problem with SOS curvature estimator for up.\n");
      }
      new_bounds = false;
      new_x = false;
      new_mults = false;

      // Determine step size and predicted change
      if (retval) {
	double alpha = 1.;
	change_up[i] = alpha*gradLagTd + 0.5*alpha*alpha*dTHLagd;
      }
      else {
	change_up[i] = -large_number;
      }

      // down (fix all members[i] with i>iWhere
      for (int j=iWhere+1; j<nMembers; j++) {
	const int& iCol = members[j];
	orig_d[iCol] = -solution[iCol];
      }
      sum = 0.;
      for (int j=0; j<=iWhere; j++) {
	sum += solution[members[j]];
      }
      d = (1.-sum)/((double)(iWhere+1));
      for (int j=0; j<=iWhere; j++) {
	orig_d[members[j]] = d;
      }

      retval = cur_estimator_->ComputeNullSpaceCurvature(
          new_bounds, numCols, solution, new_x, z_L, z_U,
	  numRows, lam, new_mults, orig_d, projected_d, gradLagTd, dTHLagd);
      // ToDo
      if (!retval) {
	printf("Problem with SOS curvature estimator for down.\n");
      }

      // Determine step size and predicted change
      if (retval) {
	double alpha = 1.;
	change_down[i] = alpha*gradLagTd + 0.5*alpha*alpha*dTHLagd;
      }
      else {
	change_down[i] = -large_number;
      }

      for (int j=0; j<nMembers; j++) {
	orig_d[members[j]] = 0.;
      }
    }
  }

  delete [] orig_d;
  delete [] projected_d;

  return -1; // We never detect an infeasible problem
}

}

// Copyright (C) 2007, International Business Machines Corporation and others. 
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Joao P. Goncalves, International Business Machines Corporation
//
// Date : November 12, 2007

#include "BonHeuristicDive.hpp"
#include "CoinHelperFunctions.hpp"
#include "CbcModel.hpp"

#include "OsiAuxInfo.hpp"

#include "CoinTime.hpp"

#include <fstream>

#include <iomanip>

using namespace std;

namespace Bonmin
{
  HeuristicDive::HeuristicDive()
    :
    CbcHeuristic(),
    setup_(NULL)
  {}

  HeuristicDive::HeuristicDive(BabSetupBase * setup)
    :
    CbcHeuristic(),
    setup_(setup)
  {
    //    Initialize(setup->options());
  }

  HeuristicDive::HeuristicDive(const HeuristicDive &copy)
    :
    CbcHeuristic(copy),
    setup_(copy.setup_)
  {}

  HeuristicDive &
  HeuristicDive::operator=(const HeuristicDive & rhs)
  {
    if(this != &rhs) {
      CbcHeuristic::operator=(rhs);
      setup_ = rhs.setup_;
    }
    return *this;
  }

  int
  HeuristicDive::solution(double &solutionValue, double *betterSolution)
  {
    if(model_->getNodeCount() || model_->getCurrentPassNumber() > 1) return 0;

    int returnCode = 0; // 0 means it didn't find a feasible solution

    OsiTMINLPInterface * nlp = dynamic_cast<OsiTMINLPInterface *>
                               (setup_->nonlinearSolver()->clone());
    TMINLP2TNLP* minlp = nlp->problem();

    // set tolerances
    double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
    double primalTolerance = 1.0e-6;

    int numberColumns;
    int numberRows;
    int nnz_jac_g;
    int nnz_h_lag;
    TNLP::IndexStyleEnum index_style;
    minlp->get_nlp_info(numberColumns, numberRows, nnz_jac_g,
			nnz_h_lag, index_style);

    const Bonmin::TMINLP::VariableType* variableType = minlp->var_types();
    const double* x_sol = minlp->x_sol();
    const double* x_l = minlp->x_l();
    const double* x_u = minlp->x_u();
    const double* g_sol = minlp->g_sol();
    const double* g_l = minlp->g_l();
    const double* g_u = minlp->g_u();

    adjustPrimalTolerance(minlp, primalTolerance);

    assert(isNlpFeasible(minlp, primalTolerance));

    // Get solution array for heuristic solution
    double* newSolution = new double [numberColumns];
    memcpy(newSolution,x_sol,numberColumns*sizeof(double));
    double* new_g_sol = new double [numberRows];


    // create a set with the indices of the fractional variables
    vector<int> integerColumns; // stores the integer variables
    int numberFractionalVariables = 0;
    for (int iColumn=0;iColumn<numberColumns;iColumn++) {
      if (variableType[iColumn] != Bonmin::TMINLP::CONTINUOUS) {
	integerColumns.push_back(iColumn);
	double value=newSolution[iColumn];
	if (fabs(floor(value+0.5)-value)>integerTolerance) {
	  numberFractionalVariables++;
	}
      }
    }

    setInternalVariables(minlp);

    int iteration = -1;
    while(numberFractionalVariables) {
      iteration++;

      // select a fractional variable to bound
      double smallestFraction = DBL_MAX;
      int bestColumn = -1;
      int bestRound = -1; // -1 rounds down, +1 rounds up
      selectVariableToBranch(minlp, integerColumns, newSolution,
			     bestColumn, bestRound);

      if(bestColumn >= 0) {
	if(bestRound < 0)
	  minlp->SetVariableUpperBound(bestColumn, floor(newSolution[bestColumn]));
	else
	  minlp->SetVariableLowerBound(bestColumn, ceil(newSolution[bestColumn]));
      } else {
	break;
      }

      nlp->initialSolve();

      if(minlp->optimization_status() != SUCCESS) {
	break;
      }

      memcpy(newSolution,x_sol,numberColumns*sizeof(double));

      numberFractionalVariables = 0;
      for(int iIntCol=0; iIntCol<(int)integerColumns.size(); iIntCol++) {
	int iColumn = integerColumns[iIntCol];
	double value=newSolution[iColumn];
	if (fabs(floor(value+0.5)-value)>integerTolerance)
	  numberFractionalVariables++;
      }

      double newSolutionValue;
      minlp->eval_f(numberColumns, newSolution, true, newSolutionValue); 
    }

    bool feasible = true;
    for (int iColumn=0;iColumn<numberColumns;iColumn++) {
      double value=newSolution[iColumn];
      if(value < x_l[iColumn] || value > x_u[iColumn]) {
	feasible = false;
	break;
      }
      if (variableType[iColumn] != Bonmin::TMINLP::CONTINUOUS) {
	if (fabs(floor(value+0.5)-value)>integerTolerance) {
	  feasible = false;
	  break;
	}
      }
    }
    minlp->eval_g(numberColumns, newSolution, true,
		  numberRows, new_g_sol);
    for(int iRow=0; iRow<numberRows; iRow++) {
      if(new_g_sol[iRow]<g_l[iRow]-primalTolerance ||
	 new_g_sol[iRow]>g_u[iRow]+primalTolerance) {
	if(minlp->optimization_status() != SUCCESS) {
	  feasible = false;
	  break;
	} else {
	  cout<<"It should be infeasible because: "<<endl;
	  cout<<"g_l["<<iRow<<"]= "<<g_l[iRow]<<" "
	      <<"g_sol["<<iRow<<"]= "<<new_g_sol[iRow]<<" "
	      <<"g_u["<<iRow<<"]= "<<g_u[iRow]<<endl;
	}
      }
    }

    if(feasible) {
      double newSolutionValue;
      minlp->eval_f(numberColumns, newSolution, true, newSolutionValue); 
      if(newSolutionValue < solutionValue) {
	memcpy(betterSolution,newSolution,numberColumns*sizeof(double));
	solutionValue = newSolutionValue;
	returnCode = 1;
      }
    }

    delete [] newSolution;
    delete [] new_g_sol;

    return returnCode;
  }


  bool
  HeuristicDive::isNlpFeasible(TMINLP2TNLP* minlp, const double primalTolerance)
  {
    int numberColumns;
    int numberRows;
    int nnz_jac_g;
    int nnz_h_lag;
    TNLP::IndexStyleEnum index_style;
    minlp->get_nlp_info(numberColumns, numberRows, nnz_jac_g,
			nnz_h_lag, index_style);

    const double* x_sol = minlp->x_sol();
    const double* x_l = minlp->x_l();
    const double* x_u = minlp->x_u();
    const double* g_sol = minlp->g_sol();
    const double* g_l = minlp->g_l();
    const double* g_u = minlp->g_u();

    // check if the problem is feasible
    for (int iColumn=0;iColumn<numberColumns;iColumn++) {
      double value=x_sol[iColumn];
      if(value < x_l[iColumn] || value > x_u[iColumn]) {
	return false;
      }
    }
    for(int iRow=0; iRow<numberRows; iRow++) {
      if(g_sol[iRow]<g_l[iRow]-primalTolerance ||
	 g_sol[iRow]>g_u[iRow]+primalTolerance) {
	return false;
      }
    }

    return true;
  }

  void
  HeuristicDive::adjustPrimalTolerance(TMINLP2TNLP* minlp, double & primalTolerance)
  {
    int numberColumns;
    int numberRows;
    int nnz_jac_g;
    int nnz_h_lag;
    TNLP::IndexStyleEnum index_style;
    minlp->get_nlp_info(numberColumns, numberRows, nnz_jac_g,
			nnz_h_lag, index_style);

    const double* g_sol = minlp->g_sol();
    const double* g_l = minlp->g_l();
    const double* g_u = minlp->g_u();

    for(int iRow=0; iRow<numberRows; iRow++) {
      if(g_sol[iRow]<g_l[iRow]-primalTolerance) {
	primalTolerance = g_l[iRow]-g_sol[iRow];
      } else if(g_sol[iRow]>g_u[iRow]+primalTolerance) {
	primalTolerance = g_sol[iRow]-g_u[iRow];
      }
    }
  }

}


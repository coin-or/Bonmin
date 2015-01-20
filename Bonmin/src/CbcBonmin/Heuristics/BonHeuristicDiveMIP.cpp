// Copyright (C) 2007, International Business Machines Corporation and others. 
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Joao P. Goncalves, International Business Machines Corporation
//
// Date : November 12, 2007

#include "BonminConfig.h"
#include "BonHeuristicDiveMIP.hpp"
#include "CoinHelperFunctions.hpp"
#include "CbcModel.hpp"
#include "BonHeuristicDive.hpp"
#include "BonSubMipSolver.hpp"
#include "BonCbcLpStrategy.hpp"

#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#endif

#include "OsiClpSolverInterface.hpp"

#include "OsiAuxInfo.hpp"

#include "CoinTime.hpp"

#include <fstream>

#include <iomanip>

#include "CoinHelperFunctions.hpp"

//#define DEBUG_BON_HEURISTIC_DIVE_MIP

using namespace std;

namespace Bonmin
{
  HeuristicDiveMIP::HeuristicDiveMIP(BonminSetup * setup)
    :
    CbcHeuristic(),
    setup_(setup),
    howOften_(100),
    mip_(NULL)
  {
    Initialize(setup);
  }

  void
  HeuristicDiveMIP::Initialize(BonminSetup * b){
    delete mip_;
    mip_ = new SubMipSolver (*b, b->prefix());
   
  }

  HeuristicDiveMIP::HeuristicDiveMIP(const HeuristicDiveMIP &copy)
    :
    CbcHeuristic(copy),
    setup_(copy.setup_),
    howOften_(copy.howOften_),
    mip_(new SubMipSolver(*copy.mip_))
  {
  }

  HeuristicDiveMIP &
  HeuristicDiveMIP::operator=(const HeuristicDiveMIP & rhs)
  {
    if(this != &rhs) {
      CbcHeuristic::operator=(rhs);
      setup_ = rhs.setup_;
      howOften_ = rhs.howOften_;
      delete mip_;
      if(rhs.mip_)
        mip_ = new SubMipSolver(*rhs.mip_);
    }
    return *this;
  }

  HeuristicDiveMIP::~HeuristicDiveMIP(){
    delete mip_;
  }

  struct MatComp{
    const int * iRow;
    const int * jCol;
  /// Destructor
    bool operator()(int i,int j){
      return (jCol[i] < jCol[j]) || (jCol[i] == jCol[j] && iRow[i] < iRow[j]);
    }
  };


  int
  HeuristicDiveMIP::solution(double &solutionValue, double *betterSolution)
  {
    if(model_->getNodeCount() || model_->getCurrentPassNumber() > 1) return 0;
    if ((model_->getNodeCount()%howOften_)!=0||model_->getCurrentPassNumber()>1)
      return 0;
 
    int returnCode = 0; // 0 means it didn't find a feasible solution

    OsiTMINLPInterface * nlp = NULL;
    if(setup_->getAlgorithm() == B_BB)
      nlp = dynamic_cast<OsiTMINLPInterface *>(model_->solver()->clone());
    else
      nlp = dynamic_cast<OsiTMINLPInterface *>(setup_->nonlinearSolver()->clone());

    TMINLP2TNLP* minlp = nlp->problem();
 
    // set tolerances
    double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
    double primalTolerance = 1.0e-6;

    int numberColumns;
    int numberRows;
    int nnz_jac_g;
    int nnz_h_lag;
    Ipopt::TNLP::IndexStyleEnum index_style;
    minlp->get_nlp_info(numberColumns, numberRows, nnz_jac_g,
			nnz_h_lag, index_style);

    const Bonmin::TMINLP::VariableType* variableType = minlp->var_types();
    const double* x_sol = minlp->x_sol();
    const double* x_l = minlp->x_l();
    const double* x_u = minlp->x_u();
    //const double* g_sol = minlp->g_sol();
    const double* g_l = minlp->g_l();
    const double* g_u = minlp->g_u();

    adjustPrimalTolerance(minlp, primalTolerance);

    assert(isNlpFeasible(minlp, primalTolerance));

    // Get information about the linear and nonlinear part of the instance
    TMINLP* tminlp = nlp->model();
    Ipopt::TNLP::LinearityType* variableLinearNonLinear = new 
      Ipopt::TNLP::LinearityType [numberColumns];
    tminlp->get_variables_linearity(numberColumns, variableLinearNonLinear);
    vector<int> linearVariable;
    vector<int> nonlinearVariable;
    for (int iColumn=0;iColumn<numberColumns;iColumn++) {
      if (variableLinearNonLinear[iColumn]==Ipopt::TNLP::LINEAR)
	linearVariable.push_back(iColumn);
      else
	nonlinearVariable.push_back(iColumn);
    }
    size_t numberLinearColumns = linearVariable.size();
    size_t numberNonlinearColumns = nonlinearVariable.size();


    // Get the indicies of the jacobian
    // This is also a way of knowing which variables are
    // used in each row
    int* indexRow = new int[nnz_jac_g];
    int* indexCol = new int[nnz_jac_g];
    minlp->eval_jac_g(numberColumns, x_sol, false,
		      numberRows, nnz_jac_g,
		      indexRow, indexCol, 0);

    vector<int> sortedIndex(nnz_jac_g);
    CoinIotaN(sortedIndex(), nnz_jac_g, 0);
    MatComp c;
    c.iRow = indexRow;
    c.jCol = indexCol;
    std::sort(sortedIndex.begin(), sortedIndex.end(), c);

    int* row = new int[nnz_jac_g];
    int* columnStart = new int[numberColumns];
    int* columnLength = new int[numberColumns];
    CoinZeroN(columnStart, numberColumns);
    CoinZeroN(columnLength, numberColumns);
    vector<vector<int> > column(numberRows); // stores the index of
    // the variables in
    // each row
    vector<vector<int> > columnInt(numberRows); // stores the index of
    // the integer variables in
    // each row
    std::vector<int> numberColumnsLinear(numberRows, 0); // stores the number
    // of the linear variables in
    // each row
    int indexCorrection = (index_style == Ipopt::TNLP::C_STYLE) ? 0 : 1;
    int iniCol = -1;
    for(int i=0; i<nnz_jac_g; i++) {
      int thisIndexCol = indexCol[sortedIndex[i]]-indexCorrection;
      if(indexCol[sortedIndex[i]] != iniCol) {
	iniCol = indexCol[sortedIndex[i]];
	columnStart[thisIndexCol] = i;
	columnLength[thisIndexCol] = 1;
      }
      else {
	columnLength[thisIndexCol]++;
      }
      row[i] = indexRow[sortedIndex[i]]-indexCorrection;
      column[row[i]].push_back(thisIndexCol);
      if (variableType[thisIndexCol] != Bonmin::TMINLP::CONTINUOUS)
	columnInt[row[i]].push_back(thisIndexCol);
      if(variableLinearNonLinear[thisIndexCol] == Ipopt::TNLP::LINEAR)
	numberColumnsLinear[row[i]]++;
    }

    // Get solution array for heuristic solution
    double* newSolution = new double [numberColumns];
    memcpy(newSolution,x_sol,numberColumns*sizeof(double));
    double* new_g_sol = new double [numberRows];


    // create a set with the indices of the fractional variables
    vector<int> integerNonlinearColumns; // stores the integer variables
    int numberFractionalNonlinearVariables = 0;
    for (size_t iNLCol=0;iNLCol<numberNonlinearColumns;iNLCol++) {
      int iColumn = nonlinearVariable[iNLCol];
      if (variableType[iColumn] != Bonmin::TMINLP::CONTINUOUS) {
	integerNonlinearColumns.push_back(iColumn);
	double value=newSolution[iColumn];
	if (fabs(floor(value+0.5)-value)>integerTolerance) {
	  numberFractionalNonlinearVariables++;
	}
      }
    }

    setInternalVariables(minlp);

    int iteration = -1;
    while(numberFractionalNonlinearVariables) {
      iteration++;

      // select a fractional variable to bound
      int bestColumn = -1;
      int bestRound = -1; // -1 rounds down, +1 rounds up
      selectVariableToBranch(minlp, integerNonlinearColumns, newSolution,
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

      if(minlp->optimization_status() != Ipopt::SUCCESS) {
	break;
      }

      memcpy(newSolution,x_sol,numberColumns*sizeof(double));

      numberFractionalNonlinearVariables = 0;
      for(int iIntCol=0; iIntCol<(int)integerNonlinearColumns.size(); iIntCol++) {
	int iColumn = integerNonlinearColumns[iIntCol];
	double value=newSolution[iColumn];
	if (fabs(floor(value+0.5)-value)>integerTolerance)
	  numberFractionalNonlinearVariables++;
      }

      double newSolutionValue;
      minlp->eval_f(numberColumns, newSolution, true, newSolutionValue); 
    }


    // now we are going to solve a MIP with the linear part of the problem
    int numberFractionalLinearVariables = 0;
    for (size_t iLCol=0;iLCol<numberLinearColumns;iLCol++) {
      int iColumn = linearVariable[iLCol];
      if (variableType[iColumn] != Bonmin::TMINLP::CONTINUOUS) {
	double value=newSolution[iColumn];
	if (fabs(floor(value+0.5)-value)>integerTolerance) {
	  numberFractionalLinearVariables++;
	}
      }
    }

    bool feasible = true;
    if(numberFractionalLinearVariables) {
      int numberMIPRows = 0;
      int* mapRows = new int[numberRows];
      for(int iRow=0; iRow<numberRows; iRow++) {
	mapRows[iRow] = -1; // this means that there are no linear columns in this row
	if(numberColumnsLinear[iRow] > 0) {
	  mapRows[iRow] = numberMIPRows++;
	}
      }

      // set all linear variables to zero in order to compute the
      // impact of the nonlinear variables in each row
      int numberIntegerLinearColumns = 0;
      for (size_t iLCol=0;iLCol<numberLinearColumns;iLCol++) {
	int iColumn = linearVariable[iLCol];
	newSolution[iColumn] = 0.0;
	if (variableType[iColumn] != Bonmin::TMINLP::CONTINUOUS)
	  numberIntegerLinearColumns++;
      }

      double* gradient_f = new double[numberColumns];
      minlp->eval_grad_f(numberColumns,newSolution,true,gradient_f);
      // create row lower and upper bounds for MILP
      minlp->eval_g(numberColumns, newSolution, true,
		    numberRows, new_g_sol);
      double* row_lb = new double[numberMIPRows];
      double* row_ub = new double[numberMIPRows];
      for(int iRow=0; iRow<numberRows; iRow++) {
	if(mapRows[iRow] > -1) {
	  assert(mapRows[iRow] < numberMIPRows);
	  if(g_l[iRow] == (-1.0) * nlp->getInfinity())
	    row_lb[mapRows[iRow]] = g_l[iRow];
	  else
	    row_lb[mapRows[iRow]] = g_l[iRow] - new_g_sol[iRow];
	  if(g_u[iRow] == nlp->getInfinity())
	    row_ub[mapRows[iRow]] = g_u[iRow];
	  else
	    row_ub[mapRows[iRow]] = g_u[iRow] - new_g_sol[iRow];
	}
      }

      // get the jacobian so that we know the coefficients of the MILP matrix
      double* jac_g = new double [nnz_jac_g];
      minlp->eval_jac_g(numberColumns, x_sol, false,
			numberRows, nnz_jac_g,
			0, 0, jac_g);


      // Define the constraint matrix for MILP
      CoinPackedMatrix* matrix = new CoinPackedMatrix(true,0,0);
      matrix->setDimensions(numberMIPRows,0);

      // create objective function and columns lower and upper bounds for MILP
      // and create columns for matrix in MILP
      double* objective = new double[numberLinearColumns];
      double* col_lb = new double[numberLinearColumns];
      double* col_ub = new double[numberLinearColumns];
      int* indexIntegerColumn = new int[numberIntegerLinearColumns];
      int numberIndexIntegerColumn = 0;
      for (size_t iLCol=0;iLCol<numberLinearColumns;iLCol++) {
	int iColumn = linearVariable[iLCol];
	objective[iLCol] = gradient_f[iColumn];
	col_lb[iLCol] = x_l[iColumn];
	col_ub[iLCol] = x_u[iColumn];
	CoinPackedVector newRow;
        int end = columnStart[iColumn]+columnLength[iColumn];
	for (int j=columnStart[iColumn];
	     j< end;j++) {
	  int iRow = row[j];
	  newRow.insert(mapRows[iRow], jac_g[sortedIndex[j]]);
	}
	matrix->appendCol(newRow);
	if (variableType[iColumn] != Bonmin::TMINLP::CONTINUOUS)
	  indexIntegerColumn[numberIndexIntegerColumn++] = static_cast<int>(iLCol);
      }

      // load the problem to OSI
      OsiSolverInterface *si = mip_->solver();
      bool delete_si = false;
      if(si == NULL){
        si = new OsiClpSolverInterface;
        mip_->setLpSolver(si);
        delete_si = true;
      }
      CoinMessageHandler * handler = model_->messageHandler()->clone();
      si->passInMessageHandler(handler);
      si->messageHandler()->setLogLevel(0);

      si->loadProblem(*matrix, col_lb, col_ub, objective, row_lb, row_ub);
      si->setInteger(indexIntegerColumn, numberIndexIntegerColumn);
      
      mip_->optimize(DBL_MAX, 0, 60);

      if(mip_->getLastSolution()) {
	const double* solution = mip_->getLastSolution();
	assert(si->getNumCols() == static_cast<int>(numberLinearColumns));
	for (size_t iLCol=0;iLCol<numberLinearColumns;iLCol++) {
	  int iColumn = linearVariable[iLCol];
	  newSolution[iColumn] = solution[iLCol];
	}
      }
      else
	feasible = false;

      delete [] mapRows;
      delete [] row_lb;
      delete [] row_ub;
      delete [] jac_g;
      delete [] gradient_f;
      delete matrix;
      delete [] objective;
      delete [] col_lb;
      delete [] col_ub;
      delete [] indexIntegerColumn;
      if(delete_si){
        delete si;
      }
      delete handler;
    }

#if 0
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
#ifdef DEBUG_BON_HEURISTIC_DIVE_MIP
	  cout<<"It should be infeasible because: "<<endl;
	  cout<<"g_l["<<iRow<<"]= "<<g_l[iRow]<<" "
	      <<"g_sol["<<iRow<<"]= "<<new_g_sol[iRow]<<" "
	      <<"g_u["<<iRow<<"]= "<<g_u[iRow]<<endl;
	  cout<<"primalTolerance= "<<primalTolerance<<endl;
#endif
	  feasible = false;
	  break;
	}
      }
    }
#else
    if(feasible) {
      // fix the integer variables and solve the NLP
      for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	if (variableType[iColumn] != Bonmin::TMINLP::CONTINUOUS) {
	  double value=newSolution[iColumn];
	  if (fabs(floor(value+0.5)-value)>integerTolerance) {
#ifdef DEBUG_BON_HEURISTIC_DIVE_MIP
	    cout<<"It should be infeasible because: "<<endl;
	    cout<<"variable "<<iColumn<<" is not integer"<<endl;
#endif
	    feasible = false;
	    break;
	  }
	  else {
	    value=floor(newSolution[iColumn]+0.5);
	    minlp->SetVariableUpperBound(iColumn, value);
	    minlp->SetVariableLowerBound(iColumn, value);
	  }
	}
      }
      if(feasible) {
	nlp->initialSolve();
	if(minlp->optimization_status() != Ipopt::SUCCESS) {
	  feasible = false;
	}
	memcpy(newSolution,x_sol,numberColumns*sizeof(double));
      }
    }
#endif

    if(feasible) {
      double newSolutionValue;
      minlp->eval_f(numberColumns, newSolution, true, newSolutionValue); 
      if(newSolutionValue < solutionValue) {
	memcpy(betterSolution,newSolution,numberColumns*sizeof(double));
	solutionValue = newSolutionValue;
	returnCode = 1;
      }
    }

    delete [] variableLinearNonLinear;
    delete [] indexRow;
    delete [] indexCol;
    delete [] row;
    delete [] columnStart;
    delete [] columnLength;
    delete [] newSolution;
    delete [] new_g_sol;
    delete nlp;

#ifdef DEBUG_BON_HEURISTIC_DIVE_MIP
    std::cout<<"DiveMIP returnCode = "<<returnCode<<std::endl;
#endif

    return returnCode;
  }
}

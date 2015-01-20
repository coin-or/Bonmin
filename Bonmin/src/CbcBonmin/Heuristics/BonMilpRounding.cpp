// Copyright (C) 2007, International Business Machines Corporation and others. 
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami CNRS
//
// Date : May, 26 2010

#include "BonMilpRounding.hpp"
#include "CoinHelperFunctions.hpp"
#include "CbcModel.hpp"
#include "BonSubMipSolver.hpp"

#include "CoinTime.hpp"

#include <fstream>

#include <iomanip>

#include "CoinHelperFunctions.hpp"
#include "OsiClpSolverInterface.hpp"

//#define DEBUG_BON_HEURISTIC_DIVE_MIP

using namespace std;

namespace Bonmin
{
  MilpRounding::MilpRounding(BonminSetup * setup)
    :
    CbcHeuristic(),
    setup_(setup),
    howOften_(20),
    mip_(NULL)
  {
    Initialize(setup);
  }

  void
  MilpRounding::Initialize(BonminSetup * b){
    delete mip_;
    mip_ = new SubMipSolver (*b, b->prefix());
   
  }

  MilpRounding::MilpRounding(const MilpRounding &copy)
    :
    CbcHeuristic(copy),
    setup_(copy.setup_),
    howOften_(copy.howOften_),
    mip_(new SubMipSolver(*copy.mip_))
  {
  }

  MilpRounding &
  MilpRounding::operator=(const MilpRounding & rhs)
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

  MilpRounding::~MilpRounding(){
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
  MilpRounding::solution(double &solutionValue, double *betterSolution)
  {
    if(model_->getCurrentPassNumber() > 1) return 0;
    if (model_->currentDepth() > 2 && (model_->getNodeCount()%howOften_)!=0)
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
    //double primalTolerance = 1.0e-6;

    int n;
    int m;
    int nnz_jac_g;
    int nnz_h_lag;
    Ipopt::TNLP::IndexStyleEnum index_style;
    minlp->get_nlp_info(n, m, nnz_jac_g,
			nnz_h_lag, index_style);

    const Bonmin::TMINLP::VariableType* variableType = minlp->var_types();
    const double* x_sol = minlp->x_sol();
    const double* g_l = minlp->g_l();
    const double* g_u = minlp->g_u();

    const double * colsol = model_->solver()->getColSolution();


    // Get information about the linear and nonlinear part of the instance
    TMINLP* tminlp = nlp->model();
    vector<Ipopt::TNLP::LinearityType> c_lin(m);
    tminlp->get_constraints_linearity(m, c_lin());
    vector<int> c_idx(m);
    int n_lin = 0;
    for (int i=0;i<m;i++) {
      if (c_lin[i]==Ipopt::TNLP::LINEAR)
	c_idx[i] = n_lin++;
      else
	c_idx[i] = -1;
    }


    // Get the structure of the jacobian
    vector<int> indexRow(nnz_jac_g);
    vector<int> indexCol(nnz_jac_g);
    minlp->eval_jac_g(n, x_sol, false,
		      m, nnz_jac_g,
		      indexRow(), indexCol(), 0);

    // get the jacobian values 
    vector<double> jac_g(nnz_jac_g);
    minlp->eval_jac_g(n, x_sol, false,
                      m, nnz_jac_g,
                      NULL, NULL, jac_g());

    // Sort the matrix to column ordered
    vector<int> sortedIndex(nnz_jac_g);
    CoinIotaN(sortedIndex(), nnz_jac_g, 0);
    MatComp c;
    c.iRow = indexRow();
    c.jCol = indexCol();
    std::sort(sortedIndex.begin(), sortedIndex.end(), c);

    vector<int> row (nnz_jac_g);
    vector<double> value (nnz_jac_g);
    vector<int> columnStart(n,0); 
    vector<int> columnLength(n,0);
    int indexCorrection = (index_style == Ipopt::TNLP::C_STYLE) ? 0 : 1;
    int iniCol = -1;
    int nnz = 0;
    for(int i=0; i<nnz_jac_g; i++) {
      int thisIndexCol = indexCol[sortedIndex[i]]-indexCorrection;
      int thisIndexRow = c_idx[indexRow[sortedIndex[i]] - indexCorrection];
      if(thisIndexCol != iniCol) {
	iniCol = thisIndexCol;
	columnStart[thisIndexCol] = nnz;
	columnLength[thisIndexCol] = 0;
      }
      if(thisIndexRow == -1) continue;
      columnLength[thisIndexCol]++;
      row[nnz] = thisIndexRow;
      value[nnz] = jac_g[i];
      nnz++;
    }

    // Build the row lower and upper bounds
    vector<double> newRowLower(n_lin);
    vector<double> newRowUpper(n_lin);
    for(int i = 0 ; i < m ; i++){
      if(c_idx[i] == -1) continue;
      newRowLower[c_idx[i]] = g_l[i];
      newRowUpper[c_idx[i]] = g_u[i];
    }

    // Get solution array for heuristic solution
    vector<double> newSolution(n);
    std::copy(x_sol, x_sol + n, newSolution.begin());

    // Define the constraint matrix for MILP
    CoinPackedMatrix matrix(true,n_lin,n, nnz, value(), row(), columnStart(), columnLength());

      // create objective function and columns lower and upper bounds for MILP
      // and create columns for matrix in MILP
      //double alpha = 0;
      double beta = 1;
      vector<double> objective(n);
      vector<int> idxIntegers;
      idxIntegers.reserve(n);
      for(int i = 0 ; i < n ; i++){
         if(variableType[i] != Bonmin::TMINLP::CONTINUOUS){
            idxIntegers.push_back(i);
            objective[i] = beta*(1 - 2*colsol[i]);
         }
      }

#if 0
      // Get dual multipliers and build gradient of the lagrangean
      const double * duals = nlp->getRowPrice() + 2 *n;
      vector<double> grad(n, 0); 
      vector<int> indices(n, 0);
      tminlp->eval_grad_f(n, x_sol, false, grad());
      for(int i = 0 ; i < m ; i++){
        if(c_lin[i] == Ipopt::TNLP::LINEAR) continue;
        int nnz;
        tminlp->eval_grad_gi(n, x_sol, false, i, nnz, indices(), NULL);  
        tminlp->eval_grad_gi(n, x_sol, false, i, nnz, NULL, grad());
        for(int k = 0 ; k < nnz ; k++){
          objective[indices[k]] += alpha *duals[i] * grad[k];
        } 
      }
      for(int i = 0 ; i < n ; i++){
         if(variableType[i] != Bonmin::TMINLP::CONTINUOUS)
         objective[i] += alpha * grad[i];
         //if(fabs(objective[i]) < 1e-4) objective[i] = 0;
         else objective[i] = 0;
      }
      std::copy(objective.begin(), objective.end(), std::ostream_iterator<double>(std::cout, " "));
      std::cout<<std::endl;
#endif

      // load the problem to OSI
      OsiSolverInterface *si = mip_->solver();
      assert(si != NULL);
      CoinMessageHandler * handler = model_->messageHandler()->clone();
      si->passInMessageHandler(handler);
      si->messageHandler()->setLogLevel(1);

      si->loadProblem(matrix, model_->solver()->getColLower(), model_->solver()->getColUpper(), objective(), 
                      newRowLower(), newRowUpper());
      si->setInteger(idxIntegers(), static_cast<int>(idxIntegers.size()));
      si->applyCuts(noGoods);

      bool hasFractionnal = true;
      while(hasFractionnal){
        mip_->optimize(DBL_MAX, 0, 60);
        hasFractionnal = false;
#if 0
        bool feasible = false;
        if(mip_->getLastSolution()) {
  	const double* solution = mip_->getLastSolution();
          std::copy(solution, solution + n, newSolution.begin());
  	feasible = true;
  
        }

    if(feasible) {
      // fix the integer variables and solve the NLP
      // also add no good cut
      CoinPackedVector v;
      double lb = 1;
      for (int iColumn=0;iColumn<n;iColumn++) {
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
            if(value > 0.5){
              v.insert(iColumn, -1);
              lb -= value;
            }
	    minlp->SetVariableUpperBound(iColumn, value);
	    minlp->SetVariableLowerBound(iColumn, value);
	  }
	}
      }
      }
#endif
      }
      bool feasible = false;
      if(mip_->getLastSolution()) {
	const double* solution = mip_->getLastSolution();
        std::copy(solution, solution + n, newSolution.begin());
	feasible = true;

        delete handler;
      }
      

    if(feasible) {
      // fix the integer variables and solve the NLP
      // also add no good cut
      CoinPackedVector v;
      double lb = 1;
      for (int iColumn=0;iColumn<n;iColumn++) {
	if (variableType[iColumn] != Bonmin::TMINLP::CONTINUOUS) {
	  double value=newSolution[iColumn];
	  if (fabs(floor(value+0.5)-value)>integerTolerance) {
	    feasible = false;
	    break;
	  }
	  else {
	    value=floor(newSolution[iColumn]+0.5);
            if(value > 0.5){
              v.insert(iColumn, -1);
              lb -= value;
            }
	    minlp->SetVariableUpperBound(iColumn, value);
	    minlp->SetVariableLowerBound(iColumn, value);
	  }
	}
      }
      OsiRowCut c;
      c.setRow(v);
      c.setLb(lb);
      c.setUb(DBL_MAX);
      noGoods.insert(c);
      if(feasible) {
	nlp->initialSolve();
	if(minlp->optimization_status() != Ipopt::SUCCESS) {
	  feasible = false;
	}
	std::copy(x_sol,x_sol+n, newSolution.begin());
      }
    }
    if(feasible) {
      double newSolutionValue;
      minlp->eval_f(n, newSolution(), true, newSolutionValue); 
      if(newSolutionValue < solutionValue) {
        std::copy(newSolution.begin(), newSolution.end(), betterSolution);
	solutionValue = newSolutionValue;
	returnCode = 1;
      }
    }

    delete nlp;

#ifdef DEBUG_BON_HEURISTIC_DIVE_MIP
    std::cout<<"DiveMIP returnCode = "<<returnCode<<std::endl;
#endif

    return returnCode;
  }
  void
  MilpRounding::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
    roptions->SetRegisteringCategory("Primal Heuristics (undocumented)", RegisteredOptions::UndocumentedCategory);
   roptions->AddStringOption2(
     "MILP_rounding_heuristic",
     "if yes runs the heuristic",
     "no",
     "no", "don't run it",
     "yes", "runs the heuristic",
     "");
  }
}

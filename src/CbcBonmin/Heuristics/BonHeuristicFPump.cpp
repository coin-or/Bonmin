// Copyright (C) 2007, International Business Machines Corporation and others. 
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Joao P. Goncalves, International Business Machines Corporation
//
// Date : November 12, 2007

#include "BonHeuristicFPump.hpp"
#include "CoinHelperFunctions.hpp"
#include "CbcModel.hpp"

#include "OsiAuxInfo.hpp"

#include "CoinTime.hpp"

#include <fstream>

#include <iomanip>

using namespace std;

//#define DEBUG_BON_HEURISTIC_FPUMP

namespace Bonmin
{
  class score_sorter {
  public:
    //! Constructor
    score_sorter(const vector<double>& score):
      score_(score) {}
    
    bool operator() (const int x, const int y) const {
      return score_[x]>score_[y];
    }
    
  private:
    const vector<double>& score_;
  };


  HeuristicFPump::HeuristicFPump()
    :
    CbcHeuristic(),
    setup_(NULL),
    objective_norm_(1),
    enableAdvanced_(false)
  {}

  HeuristicFPump::HeuristicFPump(BonminSetup * setup)
    :
    CbcHeuristic(),
    setup_(setup),
    objective_norm_(1),
    enableAdvanced_(false)
  {
    Initialize(setup->options());
  }

  HeuristicFPump::HeuristicFPump(const HeuristicFPump &copy)
    :
    CbcHeuristic(copy),
    setup_(copy.setup_),
    objective_norm_(copy.objective_norm_),
    enableAdvanced_(copy.enableAdvanced_)
  {}

  HeuristicFPump &
  HeuristicFPump::operator=(const HeuristicFPump & rhs)
  {
    if(this != &rhs) {
      CbcHeuristic::operator=(rhs);
      setup_ = rhs.setup_;
      objective_norm_ = rhs.objective_norm_;
      enableAdvanced_ = rhs.enableAdvanced_;
    }
    return *this;
  }

  int
  HeuristicFPump::solution(double &solutionValue, double *betterSolution)
  {
    if(model_->getNodeCount() || model_->getCurrentPassNumber() > 1) return 0;

    bool integerSolutionAlreadyExists = false;
    if(model_->getSolutionCount()) {
      //      bestSolutionValue = model_->getObjValue();
      integerSolutionAlreadyExists = true;
      if(!enableAdvanced_)
        return 0;
      assert(solutionValue < 1.0e50);
    }

    const int maxNumberIterations = 200;
    const double toleranceObjectiveFP = 1.0e-5;

    int returnCode = 0; // 0 means it didn't find a feasible solution

    OsiTMINLPInterface * nlp = NULL;
    if(setup_->getAlgorithm() == B_BB)
      nlp = dynamic_cast<OsiTMINLPInterface *>(model_->solver()->clone());
    else
      nlp = dynamic_cast<OsiTMINLPInterface *>(setup_->nonlinearSolver()->clone());

    TMINLP2TNLP* minlp = nlp->problem();

    // set tolerances
    double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
    double primalTolerance;
#if 0
    OsiSolverInterface * solver = model_->solver();
    solver->getDblParam(OsiPrimalTolerance,primalTolerance);
#endif
    primalTolerance=1.0e-6;

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

#ifdef DEBUG_BON_HEURISTIC_FPUMP
    const double* g_sol = minlp->g_sol();
    const double* g_l = minlp->g_l();
    const double* g_u = minlp->g_u();
    // print bounds_info
    for(int i=0; i<numberColumns; i++)
      cout<<"x_l["<<i<<"]= "<<x_l[i]<<" "
	  <<"x_sol["<<i<<"]= "<<x_sol[i]<<" "
	  <<"x_u["<<i<<"]= "<<x_u[i]<<" "
	  <<"variableType = "<<variableType[i]<<endl;
    for(int i=0; i<numberRows; i++)
      cout<<"g_l["<<i<<"]= "<<g_l[i]<<" "
	  <<"g_sol["<<i<<"]= "<<g_sol[i]<<" "
	  <<"g_u["<<i<<"]= "<<g_u[i]<<endl;

    cout<<"obj_value = "<<minlp->obj_value()<<endl;
  
    cout<<"optimization_status = "<<minlp->optimization_status()<<endl;
#endif

    // exit if the current NLP solution is infeasible
    // infeasibility is determined by the NLP solver
    if(minlp->optimization_status() != Ipopt::SUCCESS){
      delete nlp;
      return returnCode;
    }

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
    int numberIntegerColumns = (int) integerColumns.size();

    // create space to store old solutions in order to prevent cycling
    const int numberOldSolutionsStored = 4;
    double ** oldSolution = new double * [numberOldSolutionsStored];
    for (int j=0;j<numberOldSolutionsStored;j++) {
      oldSolution[j]= new double[numberIntegerColumns];
      for (int i=0;i<numberIntegerColumns;i++)
	oldSolution[j][i]=-COIN_DBL_MAX;
    }

    RoundingFPump roundObj(minlp);

    //bool stopDueToAlmostZeroObjective = false;
    double* x_bar = new double[numberIntegerColumns];
    int* indexes_x_bar = new int[numberIntegerColumns];
    double* copy_newSolution = new double[numberColumns];
    int iteration = 0;
    while(numberFractionalVariables) {
      iteration++;
      if(iteration > maxNumberIterations) {
	break;
      }
      memcpy(copy_newSolution, newSolution, numberColumns*sizeof(double));
      roundObj.round(integerTolerance, primalTolerance, copy_newSolution);
      bool flip = true;
      for(int iIntCol=0; iIntCol<numberIntegerColumns; iIntCol++) {
	int iColumn = integerColumns[iIntCol];
	double value=copy_newSolution[iColumn];
#if 0
	double value=newSolution[iColumn];
	if (fabs(floor(value+0.5)-value)>integerTolerance) {
	  value = floor(value+0.5);
	  // make sure that the new value is within bounds
	  if(value < x_l[iColumn]-primalTolerance)
	    value++;
	  else if(value > x_u[iColumn]+primalTolerance)
	    value--;
	}
#endif
	x_bar[iIntCol]=value;
	indexes_x_bar[iIntCol]=iColumn;
	if(flip && 
	   fabs(x_bar[iIntCol]-oldSolution[0][iIntCol])>integerTolerance)
	  flip = false;
      }

#ifdef DEBUG_BON_HEURISTIC_FPUMP
      cout<<"iteration= "<<iteration<<", flip= "<<flip<<endl;
#endif

      // flip some of the integer variables if the rounded solution is the
      // same as the previous one
      if(flip) {
	vector<int> sortedIntegerColumns(numberIntegerColumns);
	vector<double> score(numberIntegerColumns);
	for(int iIntCol=0; iIntCol<numberIntegerColumns; iIntCol++) {
	  int iColumn = integerColumns[iIntCol];
	  sortedIntegerColumns[iIntCol] = iIntCol;
	  double value=newSolution[iColumn];
	  score[iIntCol] = fabs(value-oldSolution[0][iIntCol]);
	}
	sort(sortedIntegerColumns.begin(),sortedIntegerColumns.end(),
	     score_sorter(score));

	int maxNumberToMove = 1;
	int numberMoved = 0;
	for(int i=0; i<numberIntegerColumns; i++) {
	  int iIntCol = sortedIntegerColumns[i];
	  if(score[iIntCol] > 0.00) {
	    int iColumn = integerColumns[iIntCol];
	    double value=newSolution[iColumn];
	    if(value-oldSolution[0][iIntCol]>0.0)
	      value = oldSolution[0][iIntCol]+1.0;
	    else
	      value = oldSolution[0][iIntCol]-1.0;
	    // make sure that the new value is within bounds
	    if(value < x_l[iColumn]-primalTolerance)
	      value++;
	    else if(value > x_u[iColumn]+primalTolerance)
	      value--;
	    assert(fabs(floor(value+0.5)-value)<=integerTolerance);
	    x_bar[iIntCol]=value;
	    numberMoved++;
	  } else
	    break;
	  if(numberMoved >= maxNumberToMove)
	    break;
	}

	// check for loop.
	// see if the new rounded solution is equal to an old solution.
	// If yes, then perturb the new rounded solution
	bool matched;
	for (int k = numberOldSolutionsStored-1; k > 0; k--) {
	  double * b = oldSolution[k];
	  matched = true;
	  for(int iIntCol=0; iIntCol<numberIntegerColumns; iIntCol++) {
	    if (fabs(x_bar[iIntCol]-b[iIntCol])>integerTolerance) {
	      matched=false;
	      break;
	    } 
	  }
	  if (matched) break;
	}

#ifdef DEBUG_BON_HEURISTIC_FPUMP
	cout<<"matched= "<<matched<<endl;
#endif

	if (matched) {
	  // perturbation
	  for(int iIntCol=0; iIntCol<numberIntegerColumns; iIntCol++) {
	    int iColumn = integerColumns[iIntCol];
	    double value=newSolution[iColumn];
	    double random = max(0.0,CoinDrand48()-0.3);
	    double difference = fabs(value-oldSolution[0][iIntCol]);
	    if(difference+random>0.5) {
	      if(value-oldSolution[0][iIntCol]>0.0)
		value = oldSolution[0][iIntCol]+1.0;
	      else
		value = oldSolution[0][iIntCol]-1.0;
	      // make sure that the new value is within bounds
	      if(value < x_l[iColumn]-primalTolerance)
		value++;
	      else if(value > x_u[iColumn]+primalTolerance)
		value--;
	      assert(fabs(floor(value+0.5)-value)<=integerTolerance);
	    } else {
	      // this variable is not going to be perturbed
	      value = oldSolution[0][iIntCol];
	    }
	    x_bar[iIntCol]=value;
	  }
	}
      }
      // store the new solution and remove the oldest one
      for (int j=numberOldSolutionsStored-1;j>0;j--) {
	for (int i = 0; i < numberIntegerColumns; i++) 
	  oldSolution[j][i]=oldSolution[j-1][i];
      }
      for (int j = 0; j < numberIntegerColumns; j++) 
	oldSolution[0][j] = x_bar[j];


      // solve the NLP problem
      double obj_nlp;
      if(integerSolutionAlreadyExists)
	// use cutoff constraint
	obj_nlp = nlp->solveFeasibilityProblem(numberIntegerColumns,
					       x_bar,indexes_x_bar,
					       objective_norm_, solutionValue);
      else
	obj_nlp = nlp->solveFeasibilityProblem(numberIntegerColumns,
					       x_bar,indexes_x_bar,
					       1,0,objective_norm_);


#ifdef DEBUG_BON_HEURISTIC_FPUMP
      cout<<"obj_nlp= "<<obj_nlp<<endl;
#endif

      memcpy(newSolution,x_sol,numberColumns*sizeof(double));

      if(obj_nlp < toleranceObjectiveFP) {
	//stopDueToAlmostZeroObjective = true;
	break;
      }

      // compute number of fractional variables
      numberFractionalVariables = 0;
      for(int iIntCol=0; iIntCol<numberIntegerColumns; iIntCol++) {
	int iColumn = integerColumns[iIntCol];
	double value=newSolution[iColumn];
	if (fabs(floor(value+0.5)-value)>integerTolerance)
	  numberFractionalVariables++;
      }

    }

    for (int j=0;j<numberOldSolutionsStored;j++) 
      delete [] oldSolution[j];
    delete [] oldSolution;
    delete [] x_bar;
    delete [] indexes_x_bar;


    // fix the integer variables and solve the NLP
    for(int iIntCol=0; iIntCol<numberIntegerColumns; iIntCol++) {
      int iColumn = integerColumns[iIntCol];
      double value=floor(newSolution[iColumn]+0.5);
      minlp->SetVariableUpperBound(iColumn, floor(value));
      minlp->SetVariableLowerBound(iColumn, ceil(value));
    }
    nlp->initialSolve();
    bool feasible = true;
    if(minlp->optimization_status() != Ipopt::SUCCESS) {
      feasible = false;
      //      if(stopDueToAlmostZeroObjective)
	//	returnCode = 8;
    }
    memcpy(newSolution,x_sol,numberColumns*sizeof(double));

    if(feasible) {
      double newSolutionValue;
      minlp->eval_f(numberColumns, newSolution, true, newSolutionValue);
      if(newSolutionValue < solutionValue) {
	memcpy(betterSolution,newSolution,numberColumns*sizeof(double));
	solutionValue = newSolutionValue;
	returnCode = 1;
      }
    }

#ifdef DEBUG_BON_HEURISTIC_FPUMP
    cout<<"returnCode= "<<returnCode<<endl;
#endif

#if 0
    delete [] indexRow;
    delete [] indexCol;
    delete [] row;
    delete [] columnStart;
    delete [] columnLength;
#endif
    delete [] newSolution;
    delete [] new_g_sol;
    delete [] copy_newSolution;
    delete nlp;

    return returnCode;
  }

  void
  HeuristicFPump::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
    roptions->SetRegisteringCategory("Primal Heuristics", RegisteredOptions::BonminCategory);
    roptions->AddBoundedIntegerOption("feasibility_pump_objective_norm","Norm of feasibility pump objective function",
				      1, 2, 1,"");
    roptions->setOptionExtraInfo("feasibility_pump_objective_norm", 63);
    roptions->AddStringOption2("heuristic_feasibility_pump", "whether the heuristic feasibility pump should be used",
      "no", "no", "", "yes", "", "");
    roptions->setOptionExtraInfo("heuristic_feasibility_pump", 63);

    roptions->SetRegisteringCategory("Primal Heuristics (undocumented)", RegisteredOptions::UndocumentedCategory);
    roptions->AddStringOption2("unstable_fp","use at your own risks",
                               "no",
                               "no", "",
                               "yes", "","");
    roptions->setOptionExtraInfo("unstable_fp", 63);
  }

  void 
  HeuristicFPump::Initialize(Ipopt::SmartPtr<Ipopt::OptionsList> options){
    options->GetIntegerValue("feasibility_pump_objective_norm", objective_norm_, "bonmin.");
    options->GetEnumValue("unstable_fp", enableAdvanced_, "bonmin.");
  }

  RoundingFPump::RoundingFPump(TMINLP2TNLP* minlp)
    :
    minlp_(minlp)
  {
    gutsOfConstructor();
  }

  RoundingFPump::~RoundingFPump()
  {
    delete [] col_and_jac_g_;
  }

  void
  RoundingFPump::gutsOfConstructor()
  {

    int nnz_jac_g;
    int nnz_h_lag;
    Ipopt::TNLP::IndexStyleEnum index_style;
    minlp_->get_nlp_info(numberColumns_, numberRows_, nnz_jac_g,
			nnz_h_lag, index_style);
    
    const double* x_sol = minlp_->x_sol();

    // Get the indicies of the jacobian
    // This is also a way of knowing which variables are
    // used in each row
    int* indexRow = new int[nnz_jac_g];
    int* indexCol = new int[nnz_jac_g];
    minlp_->eval_jac_g(numberColumns_, x_sol, false,
		       numberRows_, nnz_jac_g,
		       indexRow, indexCol, 0);

    // get the jacobian for the solution with zeros
    double* jac_g = new double [nnz_jac_g];
    double* zero_sol = new double [numberColumns_];
    minlp_->get_starting_point(numberColumns_, 1, zero_sol, 0, NULL, NULL, numberRows_, 0, NULL);
    //memset(zero_sol, 0, numberColumns_ * sizeof(double));
    minlp_->eval_jac_g(numberColumns_, zero_sol, true,
		       numberRows_, nnz_jac_g,
		       0, 0, jac_g);

    col_and_jac_g_ = new vector<pair<int, int> >[numberRows_];

    int indexCorrection = (index_style == Ipopt::TNLP::C_STYLE) ? 0 : 1;
    for(int i=0; i<nnz_jac_g; i++) {
      int thisIndexRow = indexRow[i]-indexCorrection;      
      int thisIndexCol = indexCol[i]-indexCorrection;
      pair<int, int> value(thisIndexCol, static_cast<int>(jac_g[i]));
      col_and_jac_g_[thisIndexRow].push_back(value);
    }    

    delete [] indexRow;
    delete [] indexCol;
    delete [] jac_g;
    delete [] zero_sol;
  }

  void
  RoundingFPump::round(const double integerTolerance, 
		       const double primalTolerance,
		       double* solution)
  {
    const Bonmin::TMINLP::VariableType* variableType = minlp_->var_types();
    const double* x_l = minlp_->x_l();
    const double* x_u = minlp_->x_u();
    const double* g_l = minlp_->g_l();
    const double* g_u = minlp_->g_u();


    for(int iRow=0; iRow<numberRows_; iRow++) {
      if(g_l[iRow] == 1.0 && g_u[iRow] == 1.0) {
	bool sosConstraint = true;
	double weightedSum = 0.0;
	int counter = 1;
	vector<pair<int, int> > jac_g = col_and_jac_g_[iRow];
	for (unsigned int j=0; j<jac_g.size(); j++) {
	  int iColumn = jac_g[j].first;
	  if (solution[iColumn]>=1.0-integerTolerance ||
	      jac_g[j].second != 1.0 ||
	      variableType[iColumn] == Bonmin::TMINLP::CONTINUOUS) {
	    sosConstraint = false;
	    break;
	  }
	  else {
	    weightedSum += counter * solution[iColumn];
	    counter++;
	  }
	}
#ifdef DEBUG_BON_HEURISTIC_FPUMP
	if(sosConstraint) {
	  cout<<"weightedSum= "<<weightedSum
	      <<", numberColumns_= "<<numberColumns_<<endl;
	}
#endif
        
	if(sosConstraint) {
	  double fl = floor(weightedSum + 0.5); 
	  int indexColumnSelected = static_cast<int>(fl) - 1;
          if(indexColumnSelected < 0){//Looks like all variables have been fixed to 0
            continue;
          }
	  assert(indexColumnSelected < jac_g.size());
	  for (size_t j=0; j<jac_g.size(); j++) {
	    int iColumn = jac_g[j].first;
	    if((int)j == indexColumnSelected)
	      solution[iColumn] = 1.0;
	    else
	      solution[iColumn] = 0.0;
	  }
	}
      }
    }

    for(int iColumn=0; iColumn<numberColumns_; iColumn++) {
      if(variableType[iColumn] != Bonmin::TMINLP::CONTINUOUS) {
	double value=solution[iColumn];
	if (fabs(floor(value+0.5)-value)>integerTolerance) {
	  value = floor(value+0.5);
	  // make sure that the new value is within bounds
	  if(value < x_l[iColumn]-primalTolerance)
	    value++;
	  else if(value > x_u[iColumn]+primalTolerance)
	    value--;
	  solution[iColumn] = value;
	}
      }
    }
  }
}

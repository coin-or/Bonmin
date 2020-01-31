// (C) Copyright CNRS and International Business Machines Corporation
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, LIF Université de la Méditérannée-CNRS
// Joao Goncalves, International Business Machines Corporation
//
// Date : 06/18/2008

#include "BonHeuristicLocalBranching.hpp"
#include "CbcModel.hpp"
#include "OsiBranchingObject.hpp"

namespace Bonmin {

  /** Default constructor*/
  HeuristicLocalBranching::HeuristicLocalBranching():
    LocalSolverBasedHeuristic(),
    howOften_(100),
    numberSolutions_(0)
  {}
  /** Constructor with setup.*/
  HeuristicLocalBranching::HeuristicLocalBranching(BonminSetup * setup):
    LocalSolverBasedHeuristic(setup),
    howOften_(100),
    numberSolutions_(0)
  {}

  /** Copy constructor.*/
  HeuristicLocalBranching::HeuristicLocalBranching
  (const HeuristicLocalBranching &other):
    LocalSolverBasedHeuristic(other),
    howOften_(other.howOften_),
    numberSolutions_(other.numberSolutions_)
  {}

  HeuristicLocalBranching::~HeuristicLocalBranching(){
  }

  void 
  HeuristicLocalBranching::setModel(CbcModel * model)
  {
    model_=model;
    //  assert(model_->solver());
    validate();
  }

  void
  HeuristicLocalBranching::validate()
  {
    assert(setup_ != NULL);
    OsiTMINLPInterface * nlp = dynamic_cast<OsiTMINLPInterface *>
      (setup_->nonlinearSolver());
    TMINLP2TNLP* minlp = nlp->problem();
    int numberColumns;
    int numberRows;
    int nnz_jac_g;
    int nnz_h_lag;
    Ipopt::TNLP::IndexStyleEnum index_style;
    minlp->get_nlp_info(numberColumns, numberRows, nnz_jac_g,
			nnz_h_lag, index_style);
    const Bonmin::TMINLP::VariableType* variableType = minlp->var_types();
    const double* x_l = minlp->x_l();
    const double* x_u = minlp->x_u();

    for(int i=0; i<numberColumns; i++) {
      if ((variableType[i] != Bonmin::TMINLP::CONTINUOUS) &&
	  (x_l[i] != 0.0 || x_u[i] != 1.0)) {
	setWhen(0);
	return;
      }
    }
  }

  /** Runs heuristic*/
  int
  HeuristicLocalBranching::solution(double & objectiveValue,
			  double * newSolution)
  {
    //    if(!when() || model_->getNodeCount() || model_->getCurrentPassNumber() > 1) return 0;
    if (numberSolutions_>=model_->getSolutionCount())
      return 0;
    else
      numberSolutions_=model_->getSolutionCount();

    const double * bestSolution = model_->bestSolution();
    if (!bestSolution)
      return 0; // No solution found yet

    OsiTMINLPInterface * nlp = dynamic_cast<OsiTMINLPInterface *>
                               (setup_->nonlinearSolver()->clone());


    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();

    double* vals = new double[numberIntegers];
    int* inds = new int[numberIntegers];

    for (int i=0; i<numberIntegers; i++) {
      int iColumn = integerVariable[i];
      vals[i] = bestSolution[iColumn];
      inds[i] = iColumn;
    }

    double rhs_local_branching_constraint = numberIntegers / 2; //stefan: this should be the same as floor(numInt/2) since numInt and 2 are int's
    nlp->switchToFeasibilityProblem(numberIntegers, vals, inds, rhs_local_branching_constraint);

    int r_val = 0;
    r_val = doLocalSearch(nlp, newSolution, objectiveValue, model_->getCutoff());

    delete [] vals;
    delete [] inds;

    if(r_val > 0) numberSolutions_ = model_->getSolutionCount() + 1;

    return r_val;
  }

  void
  HeuristicLocalBranching::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
   roptions->SetRegisteringCategory("Primal Heuristics (undocumented)", RegisteredOptions::UndocumentedCategory);
   roptions->AddStringOption2(
     "heuristic_local_branching",
     "if yes runs the LocalBranching heuristic",
     "no",
     "no", "",
     "yes", "",
     "");
    roptions->setOptionExtraInfo("heuristic_local_branching", 63);
  }

   /** Initiaize using passed options.*/
   void 
   HeuristicLocalBranching::Initialize(Ipopt::SmartPtr<Ipopt::OptionsList> options){
   }
}/* ends bonmin namespace*/

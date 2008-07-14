// (C) Copyright CNRS and International Business Machines Corporation
// All Rights Reserved.
// This code is published under the Common Public License.
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
    LocalSolverBasedHeuristic(){
  }
  /** Constructor with setup.*/
  HeuristicLocalBranching::HeuristicLocalBranching(BabSetupBase * setup):
    LocalSolverBasedHeuristic(setup){
  }

  /** Copy constructor.*/
  HeuristicLocalBranching::HeuristicLocalBranching
             (const HeuristicLocalBranching &other):
    LocalSolverBasedHeuristic(other){
  }

  HeuristicLocalBranching::~HeuristicLocalBranching(){
  }

  /** Runs heuristic*/
  int
  HeuristicLocalBranching::solution(double & objectiveValue,
			  double * newSolution)
  {
    if(model_->getNodeCount() || model_->getCurrentPassNumber() > 1) return 0;

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

    double rhs_local_branching_constraint = floor(numberIntegers / 2);
    nlp->switchToFeasibilityProblem(numberIntegers, vals, inds, rhs_local_branching_constraint);

    int r_val = 0;
    r_val = doLocalSearch(nlp, newSolution, objectiveValue, model_->getCutoff());

    delete [] vals;
    delete [] inds;
    delete nlp;

    return r_val;
  }

  void
  HeuristicLocalBranching::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
   roptions->SetRegisteringCategory("Local search based heuristics", RegisteredOptions::BonminCategory);
   roptions->AddStringOption2(
     "heuristic_local_branching",
     "if yes runs the LocalBranching heuristic",
     "no",
     "no", "don't run it",
     "yes", "runs the heuristic",
     "");
  }

   /** Initiaize using passed options.*/
   void 
   HeuristicLocalBranching::Initialize(Ipopt::SmartPtr<Bonmin::OptionsList> options){
   }
}/* ends bonmin namespace*/

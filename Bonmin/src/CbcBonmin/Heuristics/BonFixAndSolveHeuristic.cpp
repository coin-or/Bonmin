// (C) Copyright CNRS
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, LIF Université de la Méditérannée-CNRS
//
// Date : 06/18/2008

#include "BonFixAndSolveHeuristic.hpp"
#include "CbcModel.hpp"
#include "OsiBranchingObject.hpp"

namespace Bonmin {

  /** Default constructor*/
  FixAndSolveHeuristic::FixAndSolveHeuristic():
    LocalSolverBasedHeuristic()
    {
  }
  /** Constructor with setup.*/
  FixAndSolveHeuristic::FixAndSolveHeuristic(BonminSetup * setup):
    LocalSolverBasedHeuristic(setup)
    {
  }

  /** Copy constructor.*/
  FixAndSolveHeuristic::FixAndSolveHeuristic
             (const FixAndSolveHeuristic &other):
    LocalSolverBasedHeuristic(other){
  }

  FixAndSolveHeuristic::~FixAndSolveHeuristic(){
  }

  /** Runs heuristic*/
  int
  FixAndSolveHeuristic::solution(double & objectiveValue,
                                 double * newSolution){
    //if(model_->getNodeCount() || model_->getCurrentPassNumber() > 1) return 0;
    if(model_->getSolutionCount() > 0) return 0;
    if(model_->getNodeCount() > 1000) return 0;
    if(model_->getNodeCount() % 100 != 0) return 0;
    int numberObjects = model_->numberObjects();
    OsiObject ** objects = model_->objects();

    OsiTMINLPInterface * nlp = dynamic_cast<OsiTMINLPInterface *>
                               (model_->solver());
    if(nlp == NULL){
       nlp = dynamic_cast<OsiTMINLPInterface *>
                         (setup_->nonlinearSolver()->clone());
     }
     else {
       nlp = dynamic_cast<OsiTMINLPInterface *>(nlp->clone());
     }

    OsiBranchingInformation info = model_->usefulInformation();
    info.solution_ = model_->getColSolution();

    int dummy;
    int nFixed = 0;
    for(int i = 0 ; i < numberObjects; i++){
      if(objects[i]->infeasibility(&info, dummy) == 0.){
         objects[i]->feasibleRegion(nlp, &info);
         nFixed ++;
      }
    }
    if(nFixed < numberObjects / 3) return 0;
    double cutoff = info.cutoff_; 
    int r_val = doLocalSearch(nlp, newSolution, objectiveValue, cutoff);
    delete nlp;
    return r_val;
  }

  void
  FixAndSolveHeuristic::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
   roptions->SetRegisteringCategory("Primal Heuristics (undocumented)", RegisteredOptions::UndocumentedCategory);
   roptions->AddStringOption2(
     "fix_and_solve_heuristic",
     "if yes runs a heuristic at root where fixes all variables integer in the continuous solution",
     "no",
     "no", "don't run it",
     "yes", "runs the heuristic",
     "");
   roptions->setOptionExtraInfo("fix_and_solve_heuristic", 63);
  }

   /** Initiaize using passed options.*/
   void 
   FixAndSolveHeuristic::Initialize(Ipopt::SmartPtr<Ipopt::OptionsList> options){
   }
}/* ends bonmin namespace*/

// (C) Copyright CNRS
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, LIF Université de la Méditérannée-CNRS
//
// Date : 06/18/2008

#include "BonDummyPump.hpp"
#include "CbcModel.hpp"
#include "OsiBranchingObject.hpp"

namespace Bonmin {

  /** Default constructor*/
  DummyPump::DummyPump():
    LocalSolverBasedHeuristic(){
  }
  /** Constructor with setup.*/
  DummyPump::DummyPump(BonminSetup * setup):
    LocalSolverBasedHeuristic(setup){
  }

  /** Copy constructor.*/
  DummyPump::DummyPump
             (const DummyPump &other):
    LocalSolverBasedHeuristic(other){
  }

  DummyPump::~DummyPump(){
  }

  /** Runs heuristic*/
  int
  DummyPump::solution(double & objectiveValue,
                                 double * newSolution){
    if(model_->getNodeCount() || model_->getCurrentPassNumber() > 1) return 0;
    //int numberObjects = model_->numberObjects();
    //OsiObject ** objects = model_->objects();
    OsiTMINLPInterface * nlp = dynamic_cast<OsiTMINLPInterface *>
                               (setup_->nonlinearSolver()->clone());

    OsiBranchingInformation info = model_->usefulInformation();
    info.solution_ = model_->getColSolution();
    int numcols = model_->getNumCols();
    vector<double> vals;
    vector<int> inds;

    for(int i = 0 ;i < numcols ; i++){
       if(nlp->isInteger(i)){
            vals.push_back(info.solution_[i]);
            inds.push_back(i);
       }
    }
    nlp->switchToFeasibilityProblem(inds.size(), vals(), inds(), 1., 0., 1);

    double cutoff = info.cutoff_; 
    int r_val = doLocalSearch(nlp, newSolution, objectiveValue, cutoff);
    delete nlp;
    return r_val;
  }

  void
  DummyPump::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
   roptions->SetRegisteringCategory("Primal Heuristics (undocumented)", RegisteredOptions::UndocumentedCategory);
   roptions->AddStringOption2(
     "dummy_pump_heuristic",
     "if yes runs a heuristic which looks like a dummy FP",
     "no",
     "no", "don't run it",
     "yes", "runs the heuristic",
     "");
   roptions->setOptionExtraInfo("dummy_pump_heuristic", 63);
  }

   /** Initiaize using passed options.*/
   void 
   DummyPump::Initialize(Ipopt::SmartPtr<Ipopt::OptionsList> options){
   }
}/* ends bonmin namespace*/

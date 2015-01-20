// (C) Copyright CNRS
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, LIF Université de la Méditérannée-CNRS
//
// Date : 02/18/2009

#include "BonPumpForMinlp.hpp"
#include "CbcModel.hpp"
#include "OsiBranchingObject.hpp"

namespace Bonmin {

  /** Default constructor*/
  PumpForMinlp::PumpForMinlp():
    LocalSolverBasedHeuristic(){
  }
  /** Constructor with setup.*/
  PumpForMinlp::PumpForMinlp(BonminSetup * setup):
    LocalSolverBasedHeuristic(setup){
    setupDefaults(setup->options());
  }

  /** Copy constructor.*/
  PumpForMinlp::PumpForMinlp
             (const PumpForMinlp &other):
    LocalSolverBasedHeuristic(other){
  }

  PumpForMinlp::~PumpForMinlp(){
  }

  /** Runs heuristic*/
  int
  PumpForMinlp::solution(double & objectiveValue,
                                 double * newSolution){
    if(model_->getNodeCount() || model_->getCurrentPassNumber() > 1) return 0;
    if(model_->getSolutionCount()) return 0;
    //int numberObjects = model_->numberObjects();
    //OsiObject ** objects = model_->objects();
    OsiTMINLPInterface * nlp = dynamic_cast<OsiTMINLPInterface *>
                               (setup_->nonlinearSolver()->clone());

    OsiBranchingInformation info = model_->usefulInformation();

    double cutoff = info.cutoff_; 
    int r_val = doLocalSearch(nlp, newSolution, objectiveValue, cutoff, "pump_for_minlp.");
    return r_val;
  }

   void
   PumpForMinlp::setupDefaults(Ipopt::SmartPtr<Ipopt::OptionsList> options){
     //int dummy;
     std::string prefix = "pump_for_minlp.";
     changeIfNotSet(options, prefix, "algorithm", "B-iFP");
     changeIfNotSet(options, prefix, "time_limit", 30.);
   }


  void
  PumpForMinlp::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
   roptions->SetRegisteringCategory("Primal Heuristics", RegisteredOptions::BonminCategory);
   roptions->AddStringOption2(
     "pump_for_minlp",
     "whether to run the feasibility pump heuristic for MINLP",
     "no",
     "no", "",
     "yes", "",
     "");
    roptions->setOptionExtraInfo("pump_for_minlp", 63);
  }

   /** Initiaize using passed options.*/
   void 
   PumpForMinlp::Initialize(Ipopt::SmartPtr<Ipopt::OptionsList> options){
   }
}/* ends bonmin namespace*/


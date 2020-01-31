// (C) Copyright CNRS
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, LIF Université de la Méditérannée-CNRS
//
// Date : 06/18/2008

#include "BonLocalSolverBasedHeuristic.hpp"
#include "BonCbc.hpp"

namespace Bonmin {
  LocalSolverBasedHeuristic::LocalSolverBasedHeuristic():
     CbcHeuristic(),
     setup_(NULL),
     time_limit_(60),
     max_number_nodes_(1000),
     max_number_solutions_(10){
  }
  LocalSolverBasedHeuristic::LocalSolverBasedHeuristic(BonminSetup * setup):
     CbcHeuristic(),
     setup_(setup),
     time_limit_(60),
     max_number_nodes_(1000),
     max_number_solutions_(10){
     Initialize(setup->options());
  }

  LocalSolverBasedHeuristic::LocalSolverBasedHeuristic(const LocalSolverBasedHeuristic & other):
    CbcHeuristic(other),
    setup_(other.setup_),
    time_limit_(other.time_limit_),
    max_number_nodes_(other.max_number_nodes_),
    max_number_solutions_(other.max_number_solutions_) {
  }

   LocalSolverBasedHeuristic::~LocalSolverBasedHeuristic(){
   }  

   LocalSolverBasedHeuristic &
   LocalSolverBasedHeuristic::operator=(const LocalSolverBasedHeuristic& rhs){
     if(this != &rhs){
        CbcHeuristic::operator=(rhs);
        setup_ = rhs.setup_;
     }
     return *this;
   }

   void
   LocalSolverBasedHeuristic::changeIfNotSet(Ipopt::SmartPtr<Ipopt::OptionsList> options, 
                                             std::string prefix,
                                             const std::string &option,
                                             const std::string &value){
     int dummy;
     if(!options->GetEnumValue(option,dummy,prefix))
       options->SetStringValue(prefix + option, value, true, true);
   
   }
   void
   LocalSolverBasedHeuristic::changeIfNotSet(Ipopt::SmartPtr<Ipopt::OptionsList> options, 
                                             std::string prefix,
                                             const std::string &option,
                                             const double &value){
     double dummy;
     if(!options->GetNumericValue(option,dummy,prefix))
       options->SetNumericValue(prefix + option, value, true, true);
   
   }
   void
   LocalSolverBasedHeuristic::changeIfNotSet(Ipopt::SmartPtr<Ipopt::OptionsList> options,
                                             std::string prefix,
                                             const std::string &option,
                                             const int &value){
     int dummy;
     if(!options->GetIntegerValue(option,dummy,prefix))
       options->SetIntegerValue(prefix + option, value, true, true);
   
   }

   void
   LocalSolverBasedHeuristic::setupDefaults(Ipopt::SmartPtr<Ipopt::OptionsList> options){
     std::string prefix = "local_solver.";
     changeIfNotSet(options, prefix, "algorithm", "B-QG");
     changeIfNotSet(options, prefix, "variable_selection", "most-fractional");
     changeIfNotSet(options, prefix, "time_limit", 60.);
     changeIfNotSet(options, prefix, "node_limit", 1000);
     changeIfNotSet(options, prefix, "solution_limit", 5);
   }

   int
   LocalSolverBasedHeuristic::doLocalSearch(OsiTMINLPInterface * solver,
                                            double *solution,
                                            double & solValue,
                                            double cutoff,std::string prefix) const{
      BonminSetup * mysetup = setup_->clone(*solver, prefix);
      Bab bb;
      mysetup->setDoubleParameter(BabSetupBase::Cutoff, cutoff);
      mysetup->setIntParameter(BabSetupBase::NumberStrong, 0);
      bb(mysetup); 
      int r_val = 0;
      if(bb.bestSolution()){
        CoinCopyN(bb.bestSolution(), solver->getNumCols(), solution);
        solValue = bb.bestObj();
        r_val = 1;
      }
      delete mysetup;
      return r_val;
    }

   /** Register the options common to all local search based heuristics.*/
   void
   LocalSolverBasedHeuristic::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
   }

   /** Initiaize using passed options.*/
   void 
   LocalSolverBasedHeuristic::Initialize(Ipopt::SmartPtr<Ipopt::OptionsList> options){
     /** Some fancy defaults.*/
     setupDefaults(options);
   }
} /* Ends Bonmin namespace.*/


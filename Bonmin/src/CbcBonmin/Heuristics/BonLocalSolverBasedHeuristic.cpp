// (C) Copyright CNRS
// This code is published under the Common Public License.
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


   int
   LocalSolverBasedHeuristic::doLocalSearch(OsiTMINLPInterface * solver,
                                            double *solution,
                                            double & solValue,
                                            double cutoff) const{
     int algo;
     if (!setup_->options()->GetEnumValue("algorithm",algo,"local_solver.")) {   
      printf("Setting option\n");
      setup_->options()->SetStringValue("local_solver.algorithm","B-QG",true,true);
      setup_->options()->GetEnumValue("algorithm",algo,"local_solver.");
      printf("Option set to %i\n", algo);
    }
 
      BonminSetup * mysetup = setup_->clone(*solver, "local_solver.");
      Bab bb;
      mysetup->setDoubleParameter(BabSetupBase::Cutoff, cutoff);
      mysetup->setDoubleParameter(BabSetupBase::MaxTime, time_limit_);
      mysetup->setIntParameter(BabSetupBase::MaxNodes, max_number_nodes_);
      mysetup->setIntParameter(BabSetupBase::MaxSolutions, max_number_solutions_);
      bb(mysetup); 
      if(bb.bestSolution()){
        CoinCopyN(bb.bestSolution(), solver->getNumCols(), solution);
        bb.bestObj();
        return 1;
      }
      else return 0;
    }

   /** Register the options common to all local search based heuristics.*/
   void
   LocalSolverBasedHeuristic::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
    roptions->SetRegisteringCategory("Local search based heuristics", RegisteredOptions::BonminCategory);
    roptions->AddLowerBoundedNumberOption("local_search_time_limit","Time  limit for local searches in heuristics",
                                          0, false, 60, "");
    roptions->AddLowerBoundedIntegerOption("local_search_node_limit","Node limit for local searches in heuristics",
                                           0, 1000,"");
    roptions->AddLowerBoundedIntegerOption("local_search_solution_limit", "Solution limit for local searches in heuristics",
                               0, 5,"");
   }

   /** Initiaize using passed options.*/
   void 
   LocalSolverBasedHeuristic::Initialize(Ipopt::SmartPtr<Bonmin::OptionsList> options){
    options->GetNumericValue("local_search_time_limit", time_limit_, "bonmin..");
    options->GetIntegerValue("local_search_node_limit", max_number_nodes_, "bonmin.");
    options->GetIntegerValue("local_search_solution_limit", max_number_solutions_, "bonmin.");
   }
} /* Ends Bonmin namespace.*/


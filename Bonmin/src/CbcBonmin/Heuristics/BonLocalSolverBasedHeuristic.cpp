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
     time_limit_(180),
     max_number_nodes_(500){
  }
  LocalSolverBasedHeuristic::LocalSolverBasedHeuristic(BabSetupBase * setup):
     CbcHeuristic(),
     setup_(setup),
     time_limit_(180),
     max_number_nodes_(500){
     Initialize(setup->options());
  }

  LocalSolverBasedHeuristic::LocalSolverBasedHeuristic(const LocalSolverBasedHeuristic & other):
    CbcHeuristic(other),
    setup_(other.setup_),
    time_limit_(other.time_limit_),
    max_number_nodes_(other.max_number_nodes_) {
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
      BabSetupBase * mysetup = setup_->clone(*solver);
      Bab bb;
      mysetup->setDoubleParameter(BabSetupBase::Cutoff, cutoff);
      mysetup->setDoubleParameter(BabSetupBase::MaxTime, time_limit_);
      mysetup->setIntParameter(BabSetupBase::MaxNodes, max_number_nodes_);
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
   }

   /** Initiaize using passed options.*/
   void 
   LocalSolverBasedHeuristic::Initialize(Ipopt::SmartPtr<Bonmin::OptionsList> options){
    options->GetNumericValue("time_limit", time_limit_, "local_search.");
    options->GetIntegerValue("node_limit", max_number_nodes_, "local_search.");
   }
} /* Ends Bonmin namespace.*/


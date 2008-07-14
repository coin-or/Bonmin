// (C) Copyright CNRS and International Business Machines Corporation
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, LIF Université de la Méditérannée-CNRS
// Joao Goncalves, International Business Machines Corporation
//
// Date : 06/18/2008

#ifndef BonHeuristicLocalBranching_H
#define BonHeuristicLocalBranching_H
#include "BonLocalSolverBasedHeuristic.hpp"

namespace Bonmin {
  class HeuristicLocalBranching:public LocalSolverBasedHeuristic {
    public:
     /** Default constructor*/
     HeuristicLocalBranching();
    /** Constructor with setup.*/
    HeuristicLocalBranching(BabSetupBase * setup);

     /** Copy constructor.*/
     HeuristicLocalBranching(const HeuristicLocalBranching &other);
     /** Virtual constructor.*/
     virtual CbcHeuristic * clone() const{
      return new HeuristicLocalBranching(*this);
     }

     /** Destructor*/
     virtual ~HeuristicLocalBranching();

     /** Runs heuristic*/
     int solution(double & objectiveValue,
                  double * newSolution);
   /** Register the options common to all local search based heuristics.*/
   static void registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions);

   /** Initiaize using passed options.*/
   void Initialize(Ipopt::SmartPtr<Bonmin::OptionsList> options);
  };

}/* Ends Bonmin namepace.*/
#endif


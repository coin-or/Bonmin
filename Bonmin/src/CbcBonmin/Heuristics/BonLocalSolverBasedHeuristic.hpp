// (C) Copyright CNRS
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, LIF Université de la Méditérannée-CNRS
//
// Date : 06/18/2008

#ifndef BonLocalSolverBasedHeuristic_H
#define BonLocalSolverBasedHeuristic_H
#include "BonBabSetupBase.hpp"
#include "CbcHeuristic.hpp"

namespace Bonmin {
  class LocalSolverBasedHeuristic : public CbcHeuristic {
  public:
    /** Default constructor.*/
    LocalSolverBasedHeuristic();

    /** Constructor with setup.*/
    LocalSolverBasedHeuristic(BabSetupBase * setup);

    /** Copy constructor.*/
    LocalSolverBasedHeuristic(const LocalSolverBasedHeuristic & other);

    /** Destructor.*/
    ~LocalSolverBasedHeuristic();

    /** Virtual copy constructor.*/
    virtual CbcHeuristic * clone() const = 0;
 
  /// Assignment operator 
  LocalSolverBasedHeuristic & operator=(const LocalSolverBasedHeuristic& rhs);

#if 0
  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel * model){throw -1;}
#endif
  /// Resets stuff if model changes
  virtual void resetModel(CbcModel * model){
   setModel(model);
  }

  /** Change setup used for heuristic.*/
  void setSetup(BabSetupBase * setup){
    setup_ = setup;
    Initialize(setup_->options());
  }
  /** Performs heuristic  */
  virtual int solution(double & objectiveValue,
		       double * newSolution)=0;

  /** Performs heuristic which adds cuts  */
  virtual int solution(double & objectiveValue,
		       double * newSolution,
		       OsiCuts & cs) {return 0;}


   /** Do a local search based on setup and passed solver.*/
   int doLocalSearch(OsiTMINLPInterface * solver, 
                      double *solution, 
                      double & solValue,
                      double cutoff) const;

   /** Register the options common to all local search based heuristics.*/
   static void registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions);

   /** Initiaize using passed options.*/
   void Initialize(Ipopt::SmartPtr<Bonmin::OptionsList> options);
   protected:
   /** Setup to use for local searches (will make copies).*/
   BabSetupBase * setup_; 

   private:
    /** Time limit in local search.*/
    double time_limit_;
    /** maximal number of nodes in local search.*/
    int max_number_nodes_;
  };
} /** ends namespace Bonmin.*/

#endif


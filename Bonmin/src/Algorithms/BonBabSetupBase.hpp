// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/12/2007

#ifndef BabSetupBase_H
#define BabSetupBase_H

#include <string>
#include <list>
#include "CglCutGenerator.hpp"
#include "CbcHeuristic.hpp"
#include "OsiChooseVariable.hpp"
#include "BonOsiTMINLPInterface.hpp"
namespace Bonmin{
  /** A class to have all elements necessary to setup a branch-and-bound.*/
  class BabSetupBase {
public:
    /** Type for cut generation method with its frequency and string identification. */
    struct CuttingMethod{
      int frequency;
      std::string id;
      CglCutGenerator * cgl;
      bool atSolution;
      CuttingMethod():
      atSolution(false){
      }
      CuttingMethod(const CuttingMethod & other):
        frequency(other.frequency),
        id(other.id),
        cgl(other.cgl),
        atSolution(other.atSolution)
      {}
    };
    typedef std::list<CuttingMethod> CuttingMethods;
    typedef std::list<CbcHeuristic * > HeuristicMethods;
    
    /** Default strategies for processing next node. */
    enum NodeSelectionStrategy {
      bestBound = 0 /** Best bound*/,
      DFS /** Depth First Search*/,
      BFS /** Best First Search */,
      dynamic /** Dynamic strategy, see <a href="http://www.coin-or.org/Doxygen/Cbc/class_cbc_branch_dynamic_decision.html">
      CbcBranchActual.hpp </a> for explanations.*/
      };
     
    /** Parameters represented by an integer. */
    enum IntParameter{
      BabLogLevel = 0 /** Log level of main branch-and-bound*/,
      BabLogInterval/** Display information every logIntervval nodes.*/,
      MaxFailures /** Max number of failures in a branch.*/,
      FailureBehavior /** Behavior of the algorithm in the case of a failure.*/,
      MaxInfeasible /** Max number of consecutive infeasible problem in a branch
      before fathoming.*/,
      NumberStrong /** Number of candidates for strong branching.*/,
      MinReliability /** Minimum reliability before trust pseudo-costs.*/,
      NumEcpRoundsStrong /** Number of cutting plane iterations in lp strong branching.*/,
      MaxNodes /** Global node limit.*/,
      MaxSolutions /** limit on number of integer feasible solution.*/,
      MaxIterations /** Global iteration limit. */,
      SpecialOption /** Spetial option in particular for Cbc. */,
      DisableSos /** Consider or not SOS constraints.*/,
      NumberIntParam /** Dummy end to size table*/
    };

    
    /** Parameters represented by a double.*/
    enum DoubleParameter{
      CutoffDecr = 0 /** Amount by which cutoff is incremented */,
      Cutoff /** cutoff value */,
      AllowableGap /** Stop if absolute gap is less than this. */,
      AllowableFractionGap /** Stop if relative gap is less than this.*/,
      IntTol /** Integer tolerance.*/,
      MaxTime /** Global time limit. */,
      NumberDoubleParam /** Dummy end to size table*/
    };
    
    /** Default constructor. */
    BabSetupBase();
    
    /** Construct from existing tminlp. */
    BabSetupBase(BasicSetup& b, Ipopt::SmartPtr<TMINLP> tminlp);
    /** Construct from existing application.*/
    BabSetupBase(BasicSetup& b, Ipopt::SmartPtr<TNLPSolver> app);
    /** Construct from existing TMINLP interface.*/
    BabSetupBase(const OsiTMINLPInterface& nlp);


    /** Copy constructor. */
    BabSetupBase(const BabSetupBase & other);
    
    /** virtual copy constructor. */
    virtual BabSetupBase * clone() const = 0;
    
    /** Virtual destructor. */
    virtual ~BabSetupBase();
    
    /** Initialize from existing TMINLP interface (containing the options).*/
    void initialize(const OsiTMINLPInterface& nlp);
    /** Read options and initialize setup according to them using model in TMINLP.*/
    void initialize(Ipopt::SmartPtr<TMINLP> tminlp );
    
    /** @name Methods to instantiate: Registering and retrieving options and initializing everything. */
    /** @{ */
    /** Register all the options for this algorithm instance.*/
    virtual void registerOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions);
    /** Setup the defaults options for this algorithm. */
    virtual void setBabDefaultOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions) {}
    /** Read options and initialize algorithm according to them.*/
    virtual void initialize(OsiTMINLPInterface *);
    /** Register all the options for this algorithm instance.*/
    static void registerAllOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions);

    /** Get the basic options if don't already have them.*/
    virtual void defaultBasicOptions() = 0;
    
    /** Set the value for options, output...*/
    void setBasicOptions(BasicSetup &b){
      options_ = b.options();
      roptions_ = b.roptions();
      journalist_ = b.journalist();}
    /** @} */
    
    /** @name Elements of the branch-and-bound setup.*/
    /** @{ */
    /** Pointer to the non-linear solver used.*/
    OsiTMINLPInterface * nonlinearSolver()
    {return nonlinearSolver_;}
    /** Pointer to the continuous solver to use for relaxations. */
    OsiSolverInterface * linearSolver()
    { return linearSolver_;}
    /** list of cutting planes methods to apply with their frequencies. */
    CuttingMethods& cutGenerators()
    { return cutGenerators_;}
    /** list of Heuristic methods to use. */
    HeuristicMethods& heuristics()
    { return heuristics_;}
    /** branching method to use. */
    OsiChooseVariable * branchingMethod()
    { return branchingMethod_;}
    /** Node selection strategy. */
    NodeSelectionStrategy nodeSelectionMethod()
    {return nodeSelectionMethod_;}
    /** Return value of integer parameter. */
    int getIntParameter(const IntParameter &p)
    {return intParam_[p];}
    /** Return value of double parameter.*/
    double getDoubleParameter(const DoubleParameter &p)
    {return doubleParam_[p];}
    /** @} */
    
    void setNonlinearSolver(OsiTMINLPInterface * s){
      nonlinearSolver_ = s;}
    
    /** Get the values of base parameters from the options stored.*/
    void gatherParametersValues(){
      gatherParametersValues(options_);
    }
    
    /** Get the values of the base parameters from the passed options.*/
    void gatherParametersValues(Ipopt::SmartPtr<OptionsList> options);
    /** Acces storage of Journalist for output */
    Ipopt::SmartPtr<Ipopt::Journalist> journalist(){ return journalist_;}
    
    /** Acces list of Options */
    Ipopt::SmartPtr<Ipopt::OptionsList> options(){return options_;}
    
    /** Access registered Options */
    Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions(){return roptions_;}
    
protected:
      friend class BonminAmplSetup;
    /** storage of integer parameters.*/
    int intParam_[NumberIntParam];
    /** default values for int parameters.*/
    static int defaultIntParam_[NumberIntParam];
    /** storage of double parameters. */
    double doubleParam_[NumberDoubleParam];
    /** default values for double parameters. */
    static double defaultDoubleParam_[NumberDoubleParam];
    /** Storage of the non-linear solver used.*/
    OsiTMINLPInterface * nonlinearSolver_;
    /** Storage of continuous solver.*/
    OsiSolverInterface * linearSolver_;
    /** Cut generation methods. */
    CuttingMethods cutGenerators_;
    /** Heuristic methods. */
    HeuristicMethods heuristics_;
    /** Branching method.*/
    OsiChooseVariable * branchingMethod_;
    /** Node selection method.*/
    NodeSelectionStrategy nodeSelectionMethod_;
    
    /** Storage of Journalist for output */
    Ipopt::SmartPtr<Ipopt::Journalist> journalist_;
    
    /** List of Options */
    Ipopt::SmartPtr<Ipopt::OptionsList> options_;
    
    /** Registered Options */
    Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions_;
    
  };
}/* End namespace Bonmin. */
#endif
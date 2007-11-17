// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/09/2007

#ifndef BonIpoptHeuristic_HPP
#define BonIpoptHeuristic_HPP
#include "CbcHeuristic.hpp"
#include "BonOsiTMINLPInterface.hpp"
#include "CouenneProblem.hpp"
namespace Bonmin{
  /** An heuristic to call an NlpSolver if all CouenneObjects are close to be satisfied (for other integer objects, rounding is performed, is SOS are not satisfied does not run).
  */

  const double maxNlpInf_0 = 1e-5;

  class NlpSolveHeuristic : public CbcHeuristic{
public:
    /** Default constructor.*/
    NlpSolveHeuristic();
  /** Constructor with model and Ipopt problems.*/
    NlpSolveHeuristic(CbcModel & mip, OsiSolverInterface &nlp, bool cloneNlp = false, CouenneProblem * couenne = NULL);
    /** Copy constructor.*/
    NlpSolveHeuristic(const NlpSolveHeuristic &other);
    
    /** Destructor*/
    virtual ~NlpSolveHeuristic();
    
    /** Clone.*/
    virtual CbcHeuristic * clone() const;
    
    /** Assignment operator */
    NlpSolveHeuristic & operator=(const NlpSolveHeuristic &rhs);
    
    /** Set the nlp solver.*/
    void setNlp(OsiSolverInterface &nlp, bool cloneNlp = true);
    
    /** set the couenne problem to use.*/
    void setCouenneProblem(CouenneProblem *);
    /** Does nothing. */
    virtual void resetModel(CbcModel * model){}
    /** Run heuristic, return 1 if a better solution than the one passed is found and 0 otherwise.
        \argument objectiveValue Best known solution in input and value of solution found in output
        \argument newSolution Solution found by heuristic.
      \todo Find a quicker way to get to Couenne objects, store them or something
      */
    virtual int solution( double & objectiveValue, double * newSolution);
    /** set maxNlpInf. */
    void setMaxNlpInf(double value){
      maxNlpInf_ = value;}
    /** set number of nlp's solved for each given level of the tree*/
    void setNumberSolvePerLevel(int value){
      numberSolvePerLevel_ = value;}
private:
    /** Pointer to an nlp solver interface.*/
    OsiSolverInterface * nlp_;
    /** is nlp_ cloned or just a pointer?*/
    bool hasCloned_;
    /** maximum nlp infeasibility under which try to solve problem with Ipopt.*/
    double maxNlpInf_;
    /** Number of nlp's solved for each given level of the tree*/
    int numberSolvePerLevel_;
    /** Pointer to a couenne representation of the problem. */
    CouenneProblem * couenne_;
  };
}/* Ends namespace Bonmin. */

#endif


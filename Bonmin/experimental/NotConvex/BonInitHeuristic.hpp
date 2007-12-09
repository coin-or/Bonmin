// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Andreas Waechter, International Business Machines Corporation
//
// Date : 12/07/2007

#ifndef BonInitHeuristic_HPP
#define BonInitHeuristic_HPP

#include "CbcHeuristic.hpp"
#include "CouenneProblem.hpp"

namespace Bonmin{

  /** A heuristic that stores the initial solution of the NLP.  This
   *  is computed before Cbc is started, and in this way we can tell
   *  Cbc about this.
  */

  class InitHeuristic : public CbcHeuristic{
public:
    /** Constructor with model and Ipopt problems.*/
    InitHeuristic(double objValue, const double* sol, CouenneProblem& cp);
    /** Copy constructor.*/
    InitHeuristic(const InitHeuristic &other);
    
    /** Destructor*/
    virtual ~InitHeuristic();
    
    /** Clone.*/
    virtual CbcHeuristic * clone() const;
    
    /** Assignment operator */
    InitHeuristic & operator=(const InitHeuristic &rhs);
    
    virtual void resetModel(CbcModel * model){}
    /** Run heuristic, return 1 if a better solution than the one passed is found and 0 otherwise.
        \argument objectiveValue Best known solution in input and value of solution found in output
        \argument newSolution Solution found by heuristic.
      \todo Find a quicker way to get to Couenne objects, store them or something
      */
    virtual int solution(double & objectiveValue, double * newSolution);
private:
    /** Default constructor.*/
    InitHeuristic();

    /** objective function value from initial solve */
    double objValue_;

    /** point from initial solve */
    double* sol_;

    /** Size of array sol */
    int nVars_;
  };
}/* Ends namespace Bonmin. */

#endif


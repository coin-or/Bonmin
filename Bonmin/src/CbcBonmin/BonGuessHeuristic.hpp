// (C) Copyright International Business Machines  2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Andreas Waechter          IBM       2007-09-01

#ifndef BonGuessHeuristic_HPP
#define BonGuessHeuristic_HPP
#include "BonOsiTMINLPInterface.hpp"

#include "CbcHeuristic.hpp"

// need to change this so that we use pseudo cost objects
#include "BonChooseVariable.hpp"

namespace Bonmin
{
  class  GuessHeuristic : public CbcHeuristic
  {
  public:
    /// Usefull constructor
    GuessHeuristic(CbcModel &model);
    ///Copy constructor
    GuessHeuristic( const GuessHeuristic &copy):
      CbcHeuristic(copy)
    {}
    /// Set ChooseMethod - we need this to get pseudo costs
    inline void setChooseMethod(BonChooseVariable * chooseMethod)
    {
      chooseMethod_ = chooseMethod;
    }

    /// heuristic method providing guess, based on pseudo costs
    virtual int solution(double &solutionValue, double *betterSolution);
    virtual int solution(double &solutionValue, double *betterSolution, OsiCuts & cs)
    {
      return solution(solutionValue, betterSolution);
    }
    virtual CbcHeuristic * clone()const
    {
      return new GuessHeuristic(*this);
    }
    virtual void resetModel(CbcModel*)
    {}
  private:
    /// Default constructor
    GuessHeuristic();
    
    /// Assignment operator 
    GuessHeuristic & operator=(const GuessHeuristic& rhs);

    /// Pointer for choose method
    /// TODO: BAD THIS IS STATIC BUT I DON"T KNOW WHAT TO DO
    static BonChooseVariable * chooseMethod_;
  };
}
#endif

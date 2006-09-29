// (C) Copyright Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 03/15/2006


#ifndef CbcOaStrategy_HPP
#define CbcOaStrategy_HPP

#include "CbcStrategy.hpp"
namespace Bonmin
{
  /** A class to pass on a CbcStrategy to OA sub-milp solver.
   This class allows to setup GMI, MIR, probing and cover cuts frequency.
  Number of variables to strong branch on and minimum number of branches on a variable before its
  pseudo-cost is to be trusted.*/
  class CbcOaStrategy : public CbcStrategy
  {
  public:
    /// Default constructor
    CbcOaStrategy(
      int migFreq = -5,
      int probFreq = -5,
      int mirFreq = -5,
      int coverFreq = -5,
      int minReliability = 8,
      int numberStrong = 20,
      int nodeSelection = 0,
      double intTol = 1e-05,
      int logLevel = 0);
    /// Destructor
    virtual ~CbcOaStrategy()
    {}

    /// Virtual copy constructor
    virtual CbcStrategy * clone () const;

    /// Setup cut generators
    virtual void setupCutGenerators(CbcModel & model);
    /// Setup heuristics
    virtual void setupHeuristics(CbcModel & model);
    /// Do printing stuff
    virtual void setupPrinting(CbcModel & model,int modelLogLevel);
    /// Other stuff e.g. strong branching and preprocessing
    virtual void setupOther(CbcModel & model);


  private:
    int migFreq_;
    int probFreq_;
    int mirFreq_;
    int coverFreq_;
    int minReliability_;
    int numberStrong_;
    int nodeSelection_;
    double intTol_;
    int logLevel_;
  };
}
#endif

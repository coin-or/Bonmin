// (C) Copyright International Business Machines 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  12/26/2006


#ifndef BonOaFeasibilityChecker_HPP
#define BonOaFeasibilityChecker_HPP
#include "BonOaDecBase.hpp"

namespace Bonmin
{
  /** Class to perform OA in its classical form.*/
  class OaFeasibilityChecker : public OaDecompositionBase
  {
  public:
    /// Usefull constructor
    OaFeasibilityChecker(OsiTMINLPInterface * nlp = NULL,
        OsiSolverInterface * si = NULL,
        double cbcCutoffIncrement_=1e-07,
        double cbcIntegerTolerance = 1e-05,
        bool leaveSiUnchanged = 0
                        );

    /// New usefull constructor
    OaFeasibilityChecker(BabSetupBase &b);
    /// Copy constructor
    OaFeasibilityChecker(const OaFeasibilityChecker &copy)
        :
        OaDecompositionBase(copy)
    {}
    /// Destructor
    ~OaFeasibilityChecker();

    void setStrategy(const CbcStrategy & strategy)
    {
      parameters_.setStrategy(strategy);
    }

    virtual CglCutGenerator * clone() const
    {
      return new OaFeasibilityChecker(*this);
    }
  protected:
    /// virtual method which performs the OA algorithm by modifying lp and nlp.
    virtual double performOa(OsiCuts & cs, solverManip &nlpManip, solverManip &lpManip,
        SubMipSolver * &subMip, OsiBabSolver * babInfo, double &cutoff) const;
    /// virutal method to decide if local search is performed
    virtual bool doLocalSearch() const
    {
      return 0;
    }
  };
}
#endif

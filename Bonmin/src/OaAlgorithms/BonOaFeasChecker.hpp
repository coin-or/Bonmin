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
#include "BonOaDecHelper.hpp"

namespace Bonmin
{
  /** Class to perform OA in its classical form.*/
  class OaFeasibilityChecker : public CglCutGenerator, public OaDecompositionHelper
  {
  public:
    /// Default constructor
    OaFeasibilityChecker();
    /// Usefull constructor
    OaFeasibilityChecker(OsiTMINLPInterface * nlp = NULL,
        OsiSolverInterface * si = NULL,
        double cbcCutoffIncrement_=1e-07,
        double cbcIntegerTolerance = 1e-05,
        bool leaveSiUnchanged = 0
                   );

    /// Copy constructor
    OaFeasibilityChecker(const OaFeasibilityChecker &copy)
        :
        OaDecompositionHelper(copy)
    {
    }
    /// Destructor
    ~OaFeasibilityChecker();

    /// Assign an OsiTMINLPInterface
    void assignNlpInterface(OsiTMINLPInterface * nlp);

    /// Assign an OsiTMINLPInterface
    void assignLpInterface(OsiSolverInterface * si);

    void setStrategy(const CbcStrategy & strategy)
    {
      parameters_.setStrategy(strategy);
    }
    /// cut generation method
    virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
        const CglTreeInfo info = CglTreeInfo()) const;

    virtual CglCutGenerator * clone() const
    {
      return new OaFeasibilityChecker(*this);
    }

  };
}
#endif

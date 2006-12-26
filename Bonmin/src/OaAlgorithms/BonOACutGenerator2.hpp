// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  05/26/2005


#ifndef BonOACutGenerator2_HPP
#define BonOACutGenerator2_HPP
#include "BonOaDecHelper.hpp"

namespace Bonmin
{
  /** Class to perform OA in its classical form.*/
  class OACutGenerator2 : public CglCutGenerator, public OaDecompositionHelper
  {
  public:
    /// Default constructor
    OACutGenerator2();
    /// Usefull constructor
    OACutGenerator2(OsiTMINLPInterface * nlp = NULL,
        OsiSolverInterface * si = NULL,
        CbcStrategy * strategy = NULL,
        double cbcCutoffIncrement_=1e-07,
        double cbcIntegerTolerance = 1e-05,
        bool solveAuxiliaryProblem = 1,
        bool leaveSiUnchanged = 0
                   );

    /// Copy constructor
    OACutGenerator2(const OACutGenerator2 &copy)
        :
        OaDecompositionHelper(copy)
    {
    }
    /// Destructor
    ~OACutGenerator2();

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
      return new OACutGenerator2(*this);
    }

  };
}
#endif

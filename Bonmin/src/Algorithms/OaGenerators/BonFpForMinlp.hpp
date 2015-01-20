// (C) Copyright CNRS 2008
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// P. Bonami, CNRS
//
// Date : 02/13/2009

#ifndef BonFpForMinlp_H
#define BonFpForMinlp_H
#include "BonOaDecBase.hpp"

namespace Bonmin{
  class MinlpFeasPump : public OaDecompositionBase{
   public:
    /// Constructor with basic setup
    MinlpFeasPump(BabSetupBase & b);

    /// Copy constructor
    MinlpFeasPump(const MinlpFeasPump &copy)
        :
        OaDecompositionBase(copy),
        subMip_(new SubMipSolver(*copy.subMip_)),
        passBound_(copy.passBound_)
    {}
    /// Destructor
    ~MinlpFeasPump();

    void setStrategy(const CbcStrategy & strategy)
    {
      parameters_.setStrategy(strategy);
    }

    virtual CglCutGenerator * clone() const
    {
      return new MinlpFeasPump(*this);
    }
    /** Register OA options.*/
    static void registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions);

  protected:
    /// virtual method which performs the OA algorithm by modifying lp and nlp.
    virtual double performOa(OsiCuts & cs, solverManip &lpManip,
                   BabInfo * babInfo, double &cutoff, const CglTreeInfo & info) const;
    /// virutal method to decide if local search is performed
    virtual bool doLocalSearch(BabInfo * babInfo) const;
    /** Put objective of MIP according to FP scheme. */
    void set_fp_objective(OsiSolverInterface &si, const double * colsol) const;
    
  private:
    SubMipSolver * subMip_;
    /** Wether or not to pass bound to master algorithm.*/
    int passBound_;
  };

}/* End Namespace.*/

#endif



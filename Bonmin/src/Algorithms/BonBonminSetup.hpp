// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/13/2007
#ifndef BonminSetup_H
#define BonminSetup_H
#include "BonBabSetupBase.hpp"
namespace Bonmin{
  /* Bonmin algorithm setup. */
  class BonminSetup : public BabSetupBase{
public:
    /** Default constructor. */
    BonminSetup();
    /** Create the setup from Basic setup and existing tminlp */
    BonminSetup(BasicSetup& b, Ipopt::SmartPtr<TMINLP> tminlp);
    /** Construct the setup from an existing nlp interface.*/
    BonminSetup(const OsiTMINLPInterface& nlpSi);
    /** Copy constructor. */
    BonminSetup(const BonminSetup & other);
    /** virtual copy constructor. */
    virtual BabSetupBase * clone() const {
      return new BonminSetup(*this);}

    /** @name Methods to instantiate: Registering and retrieving options and initializing everything. */
    /** @{ */
    /** Register all the options for this algorithm instance.*/
    virtual void registerOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions);
    /** Setup the defaults options for this algorithm. */
    virtual void setBabDefaultOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions) {}
    /** @} */
    /** Register all bonmin type executable options.*/
    static void registerAllOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions);
    /** Initialize, read options and create appropriate bonmin setup.*/
    void initializeBonmin(Ipopt::SmartPtr<TMINLP> tminlp);
    /** Initialize, read options and create appropriate bonmin setup.*/
    void initializeBonmin(const OsiTMINLPInterface& nlpSi);
    /** Get the basic options if don't already have them.*/
    virtual void defaultBasicOptions();
protected:
      /** Register standard MILP cut generators. */
      static void registerMilpCutGenerators(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions);
    /** Add milp cut generators according to options.*/
    void addMilpCutGenerators();
    /** Initialize an plain branch-and-bound.*/
    void initializeBBB();
    /** Initialize a branch-and-cut with some OA.*/
    void initializeBHyb();
  };
}/** end namespace Bonmin*/

#endif


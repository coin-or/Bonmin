// (C) Copyright CNRS and others 2010
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, Université de la Méditérannée
// Hassan Hijazi, Orange Labs
//
// Date : 05/22/2010

#ifndef RobotSetup_H
#define RobotSetup_H
#include "BonBonminSetup.hpp"
#include "BonAmplSetup.hpp"
namespace Bonmin
{
  /** algorithm setup. */
  class RobotSetup : public BonminSetup
  {
  public:
    /** Default constructor. */
    RobotSetup(const CoinMessageHandler * handler = NULL);
    /** Copy constructor. */
    RobotSetup(const RobotSetup & other);

    /** Copy but uses an other nlp.*/
    RobotSetup(const RobotSetup &setup,
                OsiTMINLPInterface &nlp);

    /** Copy but uses another nlp and algorithm.*/
    RobotSetup(const RobotSetup &setup,
                OsiTMINLPInterface &nlp,
                const std::string & prefix);
    /** virtual copy constructor. */
    virtual BonminSetup * clone() const
    {
      return new RobotSetup(*this);
    }
    /** Make a copy with solver replace by one passed .*/
    //    virtual BonminSetup *clone(OsiTMINLPInterface&nlp)const{
    //      return new RobotSetup(*this, nlp);
    //    }
    /** Make a copy with solver replace by one passed .*/
    RobotSetup *clone(OsiTMINLPInterface&nlp)const{
      return new RobotSetup(*this, nlp);
    }
    /** Make a copy but take options with different prefix.*/
    RobotSetup *clone(OsiTMINLPInterface &nlp, const std::string & prefix)const{
      return new RobotSetup(*this, nlp, prefix);
    }
    virtual ~RobotSetup()
    {}
    /** @name Methods to instantiate: Registering and retrieving options and initializing everything. */
    /** @{ */
    /** Register all the options for this algorithm instance.*/
    virtual void registerOptions();
    /** @} */
    /** Register all bonmin type executable options.*/
    static void registerAllOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions);
    /** Initialize, read options and create appropriate bonmin setup.*/
    void initialize(Ipopt::SmartPtr<TMINLP> tminlp, bool createContinuousSolver = true);
    /** Initialize, read options and create appropriate bonmin setup.*/
    void initialize(const OsiTMINLPInterface& nlpSi, bool createContinuousSolver = true);

    /** Ampl initialization*/
void initialize(char **& argv)
  {
    readOptionsFile();
    /* Read the model.*/
    Ipopt::SmartPtr<AmplTMINLP> model = new AmplTMINLP(ConstPtr(journalist()), roptions(), options(),
                                                argv, NULL, "bonmin", NULL);
    mayPrintDoc();
    initialize(Ipopt::GetRawPtr(model), true);
  }



  protected:
    /** Add nway objects*/
    void addNWays();
    /** Initialize a branch-and-with robot nway.*/
    void initializeRobot();
  };
}/** end namespace Bonmin*/

#endif


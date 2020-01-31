// (C) Copyright CNRS and others 2010
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, Université de la Méditérannée
// Hassan Hijazi, Orange Labs
//
// Date : 05/22/2010

#ifndef SepaSetup_H
#define SepaSetup_H
#include "BonBonminSetup.hpp"
#include "BonAmplSetup.hpp"
namespace Sepa
{
  /** algorithm setup. */
  class SepaSetup : public Bonmin::BonminSetup
  {
  public:
    /** Default constructor. */
    SepaSetup(const CoinMessageHandler * handler = NULL);
    /** Copy constructor. */
    SepaSetup(const SepaSetup & other);

    /** Copy but uses an other nlp.*/
    SepaSetup(const SepaSetup &setup,
                Bonmin::OsiTMINLPInterface &nlp);

    /** Copy but uses another nlp and algorithm.*/
    SepaSetup(const SepaSetup &setup,
                Bonmin::OsiTMINLPInterface &nlp,
                const std::string & prefix);
    /** virtual copy constructor. */
    virtual BabSetupBase * clone() const
    {
      return new SepaSetup(*this);
    }
    /** Make a copy with solver replace by one passed .*/
    //    virtual BabSetupBase *clone(OsiTMINLPInterface&nlp)const{
    //      return new SepaSetup(*this, nlp);
    //    }
    /** Make a copy with solver replace by one passed .*/
    SepaSetup *clone(Bonmin::OsiTMINLPInterface&nlp)const{
      return new SepaSetup(*this, nlp);
    }
    /** Make a copy but take options with different prefix.*/
    SepaSetup *clone(Bonmin::OsiTMINLPInterface &nlp, const std::string & prefix)const{
      return new SepaSetup(*this, nlp, prefix);
    }
    virtual ~SepaSetup()
    {}
    /** @name Methods to instantiate: Registering and retrieving options and initializing everything. */
    /** @{ */
    /** Register all the options for this algorithm instance.*/
    virtual void registerOptions();
    /** @} */
    /** Register all bonmin type executable options.*/
    static void registerAllOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions);
    /** Initialize, read options and create appropriate bonmin setup.*/
    void initialize(Ipopt::SmartPtr<Bonmin::TMINLP> tminlp, bool createContinuousSolver = true);
    /** Initialize, read options and create appropriate bonmin setup.*/
    void initialize(const Bonmin::OsiTMINLPInterface& nlpSi, bool createContinuousSolver = true);

    /** Ampl initialization*/
void initialize(char **& argv)
  {
    readOptionsFile();
    /* Read the model.*/
    Ipopt::SmartPtr<Bonmin::AmplTMINLP> model = new Bonmin::AmplTMINLP(ConstPtr(journalist()), roptions(), options(),
                                                argv, NULL, "bonmin", NULL);
    mayPrintDoc();
    initialize(Ipopt::GetRawPtr(model), true);
  }



  protected:
    /** Initialize a branch-and-cut with some OA.*/
    void initializeSepa();
  };
}/** end namespace Bonmin*/

#endif


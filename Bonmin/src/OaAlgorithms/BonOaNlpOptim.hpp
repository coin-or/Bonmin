// (C) Copyright Carnegie Mellon University 2005, 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  05/26/2005

#ifndef BonOaNlpOptim_HPP
#define BonOaNlpOptim_HPP
#include "CglCutGenerator.hpp"
#include "BonOsiTMINLPInterface.hpp"
#include "BonOAMessages.hpp"
namespace Bonmin
{
  /** Generate cuts for the nlp corresponding to continuous relaxation at a node.*/
  class OaNlpOptim : public CglCutGenerator
  {
  public:
    /// Default constructor
    OaNlpOptim(OsiTMINLPInterface * si = NULL,
        int maxDepth = 10);
    /// Copy constructor
    OaNlpOptim(const OaNlpOptim &copy)
        :
        nlp_(copy.nlp_),
        maxDepth_(copy.maxDepth_)
    {
      handler_ = new CoinMessageHandler();
      handler_ -> setLogLevel(copy.handler_->logLevel());
      messages_ = OaMessages();
    }

    ///Abstract constructor
    virtual CglCutGenerator * clone() const
    {
      return new OaNlpOptim(*this);
    }

    /** Desctructor */
    virtual ~OaNlpOptim()
    {
      delete handler_;
    }

    /// Assign an OsiTMINLPInterface
    void assignInterface(OsiTMINLPInterface * si);
    /// cut generation method
    virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
        const CglTreeInfo info) const;



    void setMaxDepth(int value)
    {
      maxDepth_ = value;
    }
    inline int getNSolve()
    {
      return nSolve_;
    }
    /**set log level */
    void setLogLevel(int value)
    {
      handler_->setLogLevel(value);
    }
  private:
    /// Pointer to the Ipopt interface
    OsiTMINLPInterface * nlp_;

    /** maximum depth at which generate cuts*/
    int maxDepth_;

    ///Number of NLP resolution done
    mutable int nSolve_;
    /** messages handler. */
    CoinMessageHandler * handler_;
    /** handler */
    CoinMessages messages_;
  };
}
#endif

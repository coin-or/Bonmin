// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  05/26/2005

#ifndef BonOACutGenerator_HPP
#define BonOACutGenerator_HPP
#include "CglCutGenerator.hpp"
#include "BonOsiTMINLPInterface.hpp"
#include "BonOAMessages.hpp"
namespace Bonmin
{
  class OACutGenerator : public CglCutGenerator
  {
  public:
    /// Default constructor
    OACutGenerator(OsiTMINLPInterface * si = NULL,
        int maxDepth = 10);
    /// Copy constructor
    OACutGenerator(const OACutGenerator &copy)
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
      return new OACutGenerator(*this);
    }

    /** Desctructor */
    virtual ~OACutGenerator()
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

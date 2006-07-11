// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  05/26/2005

#ifndef IpCbcOACutGenerator_HPP
#define IpCbcOACutGenerator_HPP
#include "CglCutGenerator.hpp"
#include "IpoptInterface.hpp"
#include "OaMessages.hpp"

class IpCbcOACutGenerator : public CglCutGenerator
{
public:
  /// Default constructor
  IpCbcOACutGenerator(IpoptInterface * si = NULL,
      int maxDepth = 10);
  /// Copy constructor
  IpCbcOACutGenerator(const IpCbcOACutGenerator &copy)
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
    return new IpCbcOACutGenerator(*this);
  }

  /** Desctructor */
  virtual ~IpCbcOACutGenerator()
  {
    delete handler_;
  }

  /// Assign an IpoptInterface
  void assignInterface(IpoptInterface * si);
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
  IpoptInterface * nlp_;

  /** maximum depth at which generate cuts*/
  int maxDepth_;

  ///Number of NLP resolution done
  mutable int nSolve_;
  /** messages handler. */
  CoinMessageHandler * handler_;
  /** handler */
  CoinMessages messages_;
};
#endif

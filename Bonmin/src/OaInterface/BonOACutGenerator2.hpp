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
#include "CglCutGenerator.hpp"
#include "BonOsiTMINLPInterface.hpp"
#include "BonOAMessages.hpp"
#include "CbcModel.hpp"

#include "CbcStrategy.hpp"

#include "CoinTime.hpp"

namespace Bonmin
{
  class OACutGenerator2 : public CglCutGenerator
  {
  public:
    typedef enum subSolver {Clp, Cbc, Cplex, Other};
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
        :nlp_(copy.nlp_),
        si_(copy.si_),
        cbcCutoffIncrement_(copy.cbcCutoffIncrement_),
        cbcIntegerTolerance_(copy.cbcIntegerTolerance_),
        localSearchNodeLimit_(copy.localSearchNodeLimit_),
        maxLocalSearchPerNode_(copy.maxLocalSearchPerNode_),
        maxLocalSearch_(copy.maxLocalSearch_),
        maxLocalSearchTime_(copy.maxLocalSearchTime_),
        nLocalSearch_(copy.nLocalSearch_),
        solveAuxiliaryProblem_(copy.solveAuxiliaryProblem_),
        handler_(NULL), messages_(copy.messages_),
        subMilpLogLevel_(copy.subMilpLogLevel_),
        leaveSiUnchanged_(copy.leaveSiUnchanged_),
        strategy_(NULL),
        timeBegin_(0.),
        logFrequency_(copy.logFrequency_)
    {
      handler_ = new CoinMessageHandler();
      handler_->setLogLevel(copy.handler_->logLevel());
      if (copy.strategy_)
        strategy_ = copy.strategy_->clone();
      timeBegin_ = CoinCpuTime();
    }
    /// Destructor
    ~OACutGenerator2();

    /// Assign an OsiTMINLPInterface
    void assignNlpInterface(OsiTMINLPInterface * nlp);

    /// Assign an OsiTMINLPInterface
    void assignLpInterface(OsiSolverInterface * si);

    void setStrategy(const CbcStrategy & strategy)
    {
      if (strategy_)
        delete strategy_;
      strategy_ = strategy.clone();
    }
    /// cut generation method
    virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
        const CglTreeInfo info = CglTreeInfo()) const;

    virtual CglCutGenerator * clone() const
    {
      return new OACutGenerator2(*this);
    }

    inline int getNSolve()
    {
      return nSolve_;
    }
    /// Set value for cutoff increment
    void setcbcCutoffIncrement (double value)
    {
      cbcCutoffIncrement_ = value;
    }
    /// Set value for integer tolerance
    void setcbcIntegerTolerance (double value)
    {
      cbcIntegerTolerance_ = value;
    }
    ///set max number of nodes for local search
    void setLocalSearchNodeLimit(int value)
    {
      localSearchNodeLimit_ = value;
      if (si_)
        setTheNodeLimit();
    }
    ///set max number of local searches per node
    void setMaxLocalSearchPerNode(int value)
    {
      maxLocalSearchPerNode_ = value;
    }
    ///set total max number of local searches
    void setMaxLocalSearch(int value)
    {
      maxLocalSearch_ = value;
    }

    void setMaxLocalSearchTime(double time)
    {
      maxLocalSearchTime_ = time;
    }
    /**set log level */
    void setLogLevel(int value)
    {
      handler_->setLogLevel(value);
    }
    /** Set log frequency.*/
    void setLogFrequency(double value)
    {
      logFrequency_ = value;
    }
    /**set log level */
    void setSubMilpLogLevel(int value)
    {
      subMilpLogLevel_ = value;
    }
  private:
    /// Set the node limit to the interface
    void setTheNodeLimit();
    /// Set the time limit for b&b
    void setTimeLimit(double time) const;
    /// Set the cutoff for b&b
    void setCutoff(double cutoff) const;
    /// Get bound on the solution value after doing partial local search
    double siBestObj(CbcModel * model=NULL) const;
    /// Pointer to the Ipopt interface
    OsiTMINLPInterface * nlp_;
    ///Number of NLP resolution done
    mutable int nSolve_;
    /// A linear solver
    mutable OsiSolverInterface * si_;
    /// cutoff min increase (has to be intialized trhough Cbc)
    double cbcCutoffIncrement_;
    /// integer tolerance (has to be the same as Cbc's)
    double cbcIntegerTolerance_;
    ///Max number of nodes for local search
    int localSearchNodeLimit_;
    ///Max number of local searches per node
    int maxLocalSearchPerNode_;
    ///Total max number of local searches
    int maxLocalSearch_;
    /// maximum time for local searches
    double maxLocalSearchTime_;
    ///number of local searches performed
    mutable int nLocalSearch_;
    ///set to 1 to solve  an auxiliary NLP when infeasible assignment is encountered
    bool solveAuxiliaryProblem_;
    /** messages handler. */
    CoinMessageHandler * handler_;
    /** handler */
    CoinMessages messages_;
    /** sub milp log level.*/
    int subMilpLogLevel_;
    /** Wether or not we should remove cuts at the end of the procedure */
    bool leaveSiUnchanged_;
    /** Strategy to apply when using Cbc as MILP sub-solver.*/
    CbcStrategy * strategy_;
    /** time of construction*/
    double timeBegin_;
    /** Frequency of log. */
    double logFrequency_;
  };
}
#endif

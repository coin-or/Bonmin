// (C) Copyright International Business Machines Corporation, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Andreas Waechter, International Business Machines Corporation
//                    based on BonFilterSolver.cpp
//
// Date : 07/09/2007

#ifndef BonBqpdSolver_H
#define BonBqpdSolver_H

#include "BonTNLPSolver.hpp"
#include "BonBranchingTQP.hpp"

namespace Bonmin
{
  class BqpdSolver : public TNLPSolver
  {
  public:
    friend class FilterSolver;

  class UnsolvedBqpdError: public TNLPSolver::UnsolvedError
    {
    public:
      UnsolvedBqpdError(int errorNum,
          Ipopt::SmartPtr<TMINLP2TNLP> model,
          const std::string &name):
          TNLPSolver::UnsolvedError(errorNum, model, name)
      {}
      virtual const std::string& errorName() const;

      virtual const std::string& solverName() const;
      virtual ~UnsolvedBqpdError()
      {}

    private:
      static std::string errorNames_[1];
      static std::string solverName_;
    };

    /** Fortran type for integer used in filter. */
    typedef ipfint fint;
    /** Fortran type for double.used in filter */
    typedef double real;

    virtual UnsolvedError*
    newUnsolvedError(int num,
        Ipopt::SmartPtr<TMINLP2TNLP> problem,
        std::string name)
    {
      return new UnsolvedBqpdError(num, problem, name);
    }

    ///Default constructor
    BqpdSolver(bool createEmpty = false);

    /// Constructor with passed journalist, roptions, options.
    BqpdSolver(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
        Ipopt::SmartPtr<Ipopt::OptionsList> options,
        Ipopt::SmartPtr<Ipopt::Journalist> journalist
              );

    ///destructor
    virtual ~BqpdSolver();

    /** Initialize the TNLPSolver (read options from params_file)
     */
    virtual bool Initialize(std::string params_file);

    /** Initialize the TNLPSolver (read options from istream is)
     */
    virtual bool Initialize(std::istream& is);

    /** @name Solve methods */
    //@{
    /// Solves a problem expresses as a TNLP
    virtual ReturnStatus OptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> & tnlp);

    /// Resolves a problem expresses as a TNLP
    virtual ReturnStatus ReOptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> & tnlp);

    /// Set the warm start in the solver
    virtual bool setWarmStart(const CoinWarmStart * warm,
        Ipopt::SmartPtr<TMINLP2TNLP> tnlp);

    /// Get the warm start form the solver
    virtual CoinWarmStart * getWarmStart(Ipopt::SmartPtr<TMINLP2TNLP> tnlp) const;

    virtual CoinWarmStart * getEmptyWarmStart() const;

    /** Check that warm start object is valid.*/
    virtual bool warmStartIsValid(const CoinWarmStart * ws) const;

    virtual void enableWarmStart()
    {//No options to be set
    }
    virtual void disableWarmStart()
    {//No options to be set
    }
    //@}

    /// Virtual copy constructor
    virtual SmartPtr<TNLPSolver> clone();

    /// Get the CpuTime of the last optimization.
    virtual double CPUTime()
    {
      return (Ipopt::IsValid(cached_)) ? cached_->cpuTime_: 0.;
    }

    /// Get the iteration count of the last optimization.
    virtual int IterationCount()
    {
      return 0;
    }

    /// turn off all output from the solver
    virtual void turnOffOutput()
    {
      if (Ipopt::IsValid(cached_)) cached_->iprint = 0;
    }
    /// turn on all output from the solver
    virtual void turnOnOutput()
    {
      if (Ipopt::IsValid(cached_)) cached_->iprint = 3;
    }

    /// Get the solver name
    virtual std::string & solverName()
    {
      return solverName_;
    }

    /// Register this solver options into passed roptions
    void registerOptions()
    {
      registerOptions(roptions_);
    }

    /** Error code (solver specific).*/
    virtual int errorCode() const
    {
      return -1;
    }
    /// Register this solver options into passed roptions
    static void registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions);
  private:
    /** @name Private function members. */
    /** @{ */
    /** Perform optimization using data structure in cache. */
    TNLPSolver::ReturnStatus callOptimizer();
    /** @} */

    /** Cached information for reoptimizing. */
  struct cachedInfo : public Ipopt::ReferencedObject
    {
      fint n;
      fint m;
      fint k;
      fint kmax;
      real* a;
      fint* la;
      real* x;
      real* bl;
      real* bu;
      real f;
      real* g;
      real* r;
      real* w;
      real* e;
      fint* ls;
      real* alp;
      fint* lp;
      fint mlp;
      fint peq;
      real* ws;
      fint* lws;
      fint m0de;
      fint ifail;
      fint info[1];
      fint iprint;
      fint nout;

      fint kk,ll,kkk,lll,mxws,mxlws;

      Ipopt::SmartPtr<BranchingTQP> tqp_;
      /** Elapsed CPU time in last optimization. */
      double cpuTime_;
      /** flag remembering if warm start information has been put into
      cache */
      bool use_warm_start_in_cache_;

      /** Constructor.*/
      cachedInfo()
          :
          a(NULL),
          la(NULL),
          x(NULL),
          bl(NULL),
          bu(NULL),
          g(NULL),
          r(NULL),
          w(NULL),
          e(NULL),
          ls(NULL),
          alp(NULL),
          lp(NULL),
          ws(NULL),
          lws(NULL),
          cpuTime_(0),
          use_warm_start_in_cache_(false)
      {}

      cachedInfo(const Ipopt::SmartPtr<BranchingTQP> &tqp,
          Ipopt::SmartPtr<Ipopt::OptionsList>& options):
          a(NULL),
          la(NULL),
          x(NULL),
          bl(NULL),
          bu(NULL),
          g(NULL),
          r(NULL),
          w(NULL),
          e(NULL),
          ls(NULL),
          alp(NULL),
          lp(NULL),
          ws(NULL),
          lws(NULL),
          tqp_(tqp),
          cpuTime_(0),
          use_warm_start_in_cache_(false)
      {
        initialize(tqp, options);
      }

      /** Fill data structures for filter with info from tnlp. */
      void initialize(const Ipopt::SmartPtr<BranchingTQP> &tqp,
          Ipopt::SmartPtr<Ipopt::OptionsList>& options);

      /** Optimize problem described by cache with filter.*/
      void optimize();

      /** Destructor. */
      ~cachedInfo()
      {
        delete [] a;
        delete [] la;
        delete [] x;
        delete [] bl;
        delete [] bu;
        delete [] g;
        delete [] r;
        delete [] w;
        delete [] e;
        delete [] ls;
        delete [] alp;
        delete [] lp;
        delete [] ws;
        delete [] lws;
      }
    };

    /** Cached information on last problem optimized for reoptimization. */
    Ipopt::SmartPtr<cachedInfo> cached_;

    //name of solver (Bqpd)
    static std::string  solverName_;
  };

}// end namespace Bonmin
#endif

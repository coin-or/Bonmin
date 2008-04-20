// (C) Copyright International Business Machines Corporation, 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 10/02/2006

#ifndef FilterSolver_H
#define FilterSolver_H

#include "BonTNLPSolver.hpp"
#include "BonFilterTypes.hpp"
#include "BonBqpdSolver.hpp"
#include "BonFilterWarmStart.hpp"

namespace Bonmin
{
  class FilterSolver : public TNLPSolver
  {
  public:

    friend struct BqpdSolver::cachedInfo;

  class UnsolvedFilterError: public TNLPSolver::UnsolvedError
    {
    public:
      UnsolvedFilterError(int errorNum,
          Ipopt::SmartPtr<TMINLP2TNLP> model,
          const std::string &name):
          TNLPSolver::UnsolvedError(errorNum, model, name)
      {}
      virtual const std::string& errorName() const;

      virtual const std::string& solverName() const;
      virtual ~UnsolvedFilterError()
      {}

    private:
      static std::string errorNames_[1];
      static std::string solverName_;
    };

    /** Fortran type for integer used in filter. */
    typedef FilterTypes::fint fint;
    /** Fortran type for double.used in filter */
    typedef FilterTypes::real real;


    virtual UnsolvedError * newUnsolvedError(int num,
        Ipopt::SmartPtr<TMINLP2TNLP> problem,
        std::string name)
    {
      return new UnsolvedFilterError(num, problem, name);
    }


    ///Default constructor
    FilterSolver(bool createEmpty = false);


    /// Constructor with passed journalist, roptions, options.
    FilterSolver(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
        Ipopt::SmartPtr<Ipopt::OptionsList> options,
        Ipopt::SmartPtr<Ipopt::Journalist> journalist
                );

    ///destructor
    virtual ~FilterSolver();

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

   /// Get warm start used in last optimization
   virtual CoinWarmStart * getUsedWarmStart(Ipopt::SmartPtr<TMINLP2TNLP> tnlp) const{
     if(warmF_.IsValid())
       return new FilterWarmStart(*warmF_);
     else return NULL;
   }


    /// Get the warm start form the solver
    virtual CoinWarmStart * getWarmStart(Ipopt::SmartPtr<TMINLP2TNLP> tnlp) const;

    virtual CoinWarmStart * getEmptyWarmStart() const;

    /** Check that warm start object is valid.*/
    virtual bool warmStartIsValid(const CoinWarmStart * ws) const;


    virtual void enableWarmStart()
    {//No options to be set
    }
    virtual void disableWarmStart()
    {
      warmF_ = NULL;
     //No options to be set
    }
    //@}

    /// Virtual copy constructor
    virtual SmartPtr<TNLPSolver> clone();

    /// Get the CpuTime of the last optimization.
    virtual double CPUTime()
    {
      return (cached_.IsValid()) ? cached_->cpuTime_: 0.;
    }

    /// Get the iteration count of the last optimization.
    virtual int IterationCount()
    {
      return (cached_.IsValid()) ? cached_->istat[1]:0;
    }

    /// turn off all output from the solver
    virtual void turnOffOutput()
    {
      if (cached_.IsValid()) cached_->iprint = 0;
    }
    /// turn on all output from the solver
    virtual void turnOnOutput()
    {
      if (cached_.IsValid()) cached_->iprint = 3;
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
    Coin::SmartPtr<FilterWarmStart> warmF_;

    /** Cached information for reoptimizing. */
  struct cachedInfo : public Coin::ReferencedObject
    {
      fint n;
      fint m;
      fint nnz_h_;
      fint kmax;
      fint maxa;
      fint maxf;
      fint mlp;
      fint maxWk;
      fint maxiWk;
      fint iprint;
      fint nout;
      fint ifail;
      real rho;
      real * x;
      real * c;
      real f;
      real fmin;
      real * bounds;
      real * s;
      real * a;
      fint * la;
      real * ws;
      fint * lws;
      real * lam;
      char * cstype;
      fint maxiter;
      fint * istat;
      real * rstat;
      Ipopt::TNLP * tnlp_;
      fint * hStruct_;
      int * permutationJac_;
      int * permutationHess_;
      /** Elapsed CPU time in last optimization. */
      double cpuTime_;
      /** flag remembering if warm start information has been put into
      cache */
      bool use_warm_start_in_cache_;


      /** Constructor.*/
      cachedInfo():
          n(-1),
          m(-1),
          nnz_h_(-1),
          kmax(-1),
          maxa(-1),
          maxf(-1),
          mlp(-1),
          maxWk(-1),
          maxiWk(-1),
          iprint(-1),
          nout(6),
          ifail(-100),
          rho(0),
          x(NULL),
          c(NULL),
          f(1e100),
          fmin(-1e100),
          bounds(NULL),
          s(NULL),
          a(NULL),
          la(NULL),
          ws(NULL),
          lws(NULL),
          lam(NULL),
          cstype(NULL),
          maxiter(1000),
          istat(NULL),
          rstat(NULL),
          tnlp_(NULL),
          hStruct_(NULL),
          permutationJac_(NULL),
          permutationHess_(NULL),
          cpuTime_(0),
          use_warm_start_in_cache_(false)
      {}

      cachedInfo(const Ipopt::SmartPtr<Ipopt::TNLP> &tnlp,
          Ipopt::SmartPtr<Ipopt::OptionsList>& options):
          n(-1),
          m(-1),
          nnz_h_(-1),
          kmax(-1),
          maxa(-1),
          maxf(-1),
          mlp(-1),
          maxWk(-1),
          maxiWk(-1),
          iprint(-1),
          nout(6),
          ifail(0),
          rho(0),
          x(NULL),
          c(NULL),
          f(1e100),
          fmin(-1e100),
          bounds(NULL),
          s(NULL),
          a(NULL),
          la(NULL),
          ws(NULL),
          lws(NULL),
          lam(NULL),
          cstype(NULL),
          maxiter(1000),
          istat(NULL),
          rstat(NULL),
          tnlp_(NULL),
          hStruct_(NULL),
          permutationJac_(NULL),
          permutationHess_(NULL),
          cpuTime_(0),
          use_warm_start_in_cache_(false)
      {
        initialize(tnlp, options);
      }

      /** Fill data structures for filter with info from tnlp. */
      void initialize(const Ipopt::SmartPtr<Ipopt::TNLP> &tnlp,
          Ipopt::SmartPtr<Ipopt::OptionsList>& options);

      /** Optimize problem described by cache with filter.*/
      void optimize();

      /** Destructor. */
      ~cachedInfo()
      {
        delete [] x;
        delete [] c;
        delete [] bounds;
        delete [] s;
        delete [] a;
        delete [] la;
        delete [] ws;
        delete [] lws;
        delete [] lam;
        delete [] cstype;
        delete [] istat;
        delete [] rstat;
        delete [] permutationJac_;
        delete [] permutationHess_;
        delete [] hStruct_;
        tnlp_ = NULL;
      }

      void load_ws(Coin::SmartPtr<FilterWarmStart>);
    };

    /** Cached information on last problem optimized for reoptimization. */
    Coin::SmartPtr<cachedInfo> cached_;

    //name of solver (Filter)
    static std::string  solverName_;

    /** Converting TMatrices into row-ordered matrices */
    static void TMat2RowPMat(bool symmetric, fint n, fint m, int nnz, const Index* iRow,
        const Index* iCol, int * permutation2,
        fint * lws, int nnz_offset, int n_offset,
        Ipopt::TNLP::IndexStyleEnum index_style);
  };

}// end namespace Bonmin
#endif

// (C) Copyright International Business Machines Corporation, 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
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

//#define TIME_BQPD
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
        Ipopt::SmartPtr<Ipopt::Journalist> journalist,
        const std::string & prefix
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

    /// Safe the current state (after most recent solve that must have
    /// been successful) as hot start information and use that for all
    /// further solves, until unmarkHotStart is called.
    virtual bool markHotStart(){return cached_->markHotStart();}

    /// Get warm start used in last optimization
    virtual CoinWarmStart * getUsedWarmStart(Ipopt::SmartPtr<TMINLP2TNLP> tnlp) const{
      throw CoinError(__PRETTY_FUNCTION__,"","Not implemented");
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
    {//No options to be set
    }
    //@}

    /// Virtual copy constructor
    virtual Ipopt::SmartPtr<TNLPSolver> clone();

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


    // turn off all output from the solver
    virtual void setOutputToDefault()
    {
      if (Ipopt::IsValid(cached_)) cached_->iprint = default_log_level_;
    }
    /// turn on all output from the solver
    virtual void forceSolverOutput(int log_level)
    {
      if (Ipopt::IsValid(cached_)) cached_->iprint = log_level;
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
#ifdef TIME_BQPD
  struct Times {
    int numsolve;
    double create;
    double solve;
    double resolve;
    double warm_start;
    int pivots;
    Times(): numsolve(0), create(0), solve(0), resolve(0), warm_start(0), pivots(0){
    }

    Times & operator +=(Times &rhs){
      numsolve += rhs.numsolve;
      create += rhs.create;
      solve += rhs.solve;
      resolve += rhs.resolve;
      warm_start += rhs.warm_start;
      pivots += rhs.pivots;
      return *this;
    }
  };
#endif

    /// Register this solver options into passed roptions
    static void registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions);
  private:
    /** @name Private function members. */
    /** @{ */
    /** Perform optimization using data structure in cache. */
    TNLPSolver::ReturnStatus callOptimizer();
    /** @} */

    /** @name User options */
    //@{
    /** Fill-in factor for QP factorization */
    double fillin_factor_;
    int kmax_ipt_;
    int mlp_ipt_;
    /** Hot start m0de.*/
    static int m0de_;
    /** Hot start reinitialization fequency.*/
    static int reinit_freq_;
    //@}

    /** Cached information for reoptimizing. */
  struct cachedInfo : public Ipopt::ReferencedObject
    {
#ifdef TIME_BQPD
  Times times_;
#endif

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

      /// wsc common block
      fint kk,ll,mxws,mxlws;

      /** indicates if we should start from a hotstart **/
      bool haveHotStart_;
      /** @name All remaining information from all common blocks for
	  hot start **/
      //@{
      /// bqpdc common block
      fint irh1,na,na1,nb,nb1,ka1,kb1,kc1,irg1,lu1,lv,lv1,ll1;
      /// epsc common block
      real eps,tol,emin;
      /// vstepc common block
      real vstep;
      /// repc common block
      real sgnf;
      fint nrep,npiv,nres;
      /// refactorc common block
      fint nup,nfreq;
      /// alphac common block
      real alpha;
      /// sparsec common block
      fint ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof,lc;
      fint lc1,li,li1,lm,lm1,lp_,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1;
      /// factorc common block
      fint m1,m2,mp,mq,lastr,irow;
      /// mxm1c common block
      fint mxm1;
      /// /minorc
      real c;
      //@}
      fint kHot;
      real* xHot;
      real fHot;
      real* gHot;
      real* rHot;
      real* wHot;
      real* eHot;
      fint* lsHot;
      real* alpHot;
      fint* lpHot;
      fint peqHot;
      real* wsHot;
      fint* lwsHot;
      fint infoHot[1];
      fint kkkHot;
      fint lllHot;
      fint kkHot;
      fint llHot;

      Ipopt::SmartPtr<BranchingTQP> tqp_;
      /** Elapsed CPU time in last optimization. */
      double cpuTime_;
      /** flag remembering if warm start information has been put into
      cache */
      bool use_warm_start_in_cache_;
      bool bad_warm_start_info_;

      /** Number of nonzeros in Jacobian and gradient */
      int amax_;

      /** Fill-in factor for QP factorization.  This is a pointer to
	  the corresponding value in the BqpdSolver object, so that an
	  increase is not forgotten. */
      double* fillin_factor_;
      //@}

      /** next reinitialization of hot start.*/
      int next_reinit_;
      /** Number of pivots.*/
      int pivots_;
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
	  haveHotStart_(false),
	  xHot(NULL),
	  gHot(NULL),
	  rHot(NULL),
	  wHot(NULL),
	  eHot(NULL),
	  lsHot(NULL),
	  alpHot(NULL),
	  lpHot(NULL),
	  wsHot(NULL),
	  lwsHot(NULL),
          cpuTime_(0),
          use_warm_start_in_cache_(false),
	  bad_warm_start_info_(false),
          next_reinit_(0),
          pivots_(0)
      {}

      cachedInfo(const Ipopt::SmartPtr<BranchingTQP> &tqp,
		 Ipopt::SmartPtr<Ipopt::OptionsList>& options,
		 int kmax_ipt, int mlp_ipt, double* fillin_factor):
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
	  haveHotStart_(false),
	  xHot(NULL),
	  gHot(NULL),
	  rHot(NULL),
	  wHot(NULL),
	  eHot(NULL),
	  lsHot(NULL),
	  alpHot(NULL),
	  lpHot(NULL),
	  wsHot(NULL),
	  lwsHot(NULL),
          tqp_(tqp),
          cpuTime_(0),
          use_warm_start_in_cache_(false),
          bad_warm_start_info_(false),
          next_reinit_(0),
          pivots_(0)
      {
        initialize(tqp, options, kmax_ipt, mlp_ipt, fillin_factor);
      }

      /** Fill data structures for filter with info from tnlp. */
      void initialize(const Ipopt::SmartPtr<BranchingTQP> &tqp,
		      Ipopt::SmartPtr<Ipopt::OptionsList>& options,
		      int kmax_ipt, int mlp_ipt, double* fillin_factor);

      void re_initialize();
      /** Optimize problem described by cache with filter.*/
      void optimize();

      /** Store most recent solution as hot start */
      bool markHotStart();

      /** Forget about the hot start info */
      void unmarkHotStart();

      /** Copy current values from hot start info */
      void copyFromHotStart();

      /** Destructor. */
      ~cachedInfo();
    };

    /** Cached information on last problem optimized for reoptimization. */
    Ipopt::SmartPtr<cachedInfo> cached_;

    ///\name of solver (Bqpd)
    static std::string  solverName_;

    /** To record default log level.*/
     int default_log_level_;

    public:

#ifdef TIME_BQPD
    Times &times(){ return cached_->times_;}
#endif

  };

}// end namespace Bonmin
#endif

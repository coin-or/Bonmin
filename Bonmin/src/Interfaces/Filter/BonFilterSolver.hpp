// (C) Copyright International Business Machines Corporation, Carnegie Mellon University 2006
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

#include "CoinWarmStartBasis.hpp"
namespace Bonmin{
class FilterSolver : public TNLPSolver{
public:
  class UnsolvedFilterError: public TNLPSolver::UnsolvedError
{
 public:
  UnsolvedFilterError(int errorNum = 10000):
  TNLPSolver::UnsolvedError(errorNum)
  {}
  virtual const std::string& errorName() const;
 
  virtual const std::string& solverName() const; 
  virtual ~UnsolvedFilterError(){}
  private:
  static std::string errorNames_[1];
  static std::string solverName_;
};

  virtual UnsolvedError * newUnsolvedError(int num){
    return new UnsolvedFilterError(num);}



    typedef long fint;
    typedef double real;

  ///Default constructor
  FilterSolver();

  ///destructor
  virtual ~FilterSolver();

  /** Initialize the TNLPSolver (read options from params_file)
   */
   virtual void Initialize(std::string params_file);

   /** Initialize the TNLPSolver (read options from istream is)
   */
   virtual void Initialize(std::istream& is);

   /** @name Solve methods */
   //@{
   /// Solves a problem expresses as a TNLP 
  virtual ReturnStatus OptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> & tnlp);

   /// Resolves a problem expresses as a TNLP 
   virtual ReturnStatus ReOptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> & tnlp);

  /// Set the warm start in the solver
  virtual bool setWarmStart(const CoinWarmStart * warm, 
                            Ipopt::SmartPtr<TMINLP2TNLP> tnlp){
   return 1;}

  /// Get the warm start form the solver
  virtual CoinWarmStart * getWarmStart(Ipopt::SmartPtr<TMINLP2TNLP> tnlp) const{
   return new CoinWarmStartBasis;}

  virtual CoinWarmStart * getEmptyWarmStart() const {
  return new CoinWarmStartBasis;}


  virtual void enableWarmStart()
  {if(Ipopt::IsValid(cached_)) cached_->ifail= -1;}
  virtual void disableWarmStart()
  {if(Ipopt::IsValid(cached_)) cached_->ifail= 0;}
   //@}

  /// Virtual constructor
  virtual TNLPSolver * createNew(){
    return new FilterSolver;}

  ///Get a pointer to a journalist
  virtual Ipopt::SmartPtr<Ipopt::Journalist> Jnlst(){
    return journalist_;}

  ///Get a pointer to RegisteredOptions (generally used to add new ones)
  virtual Ipopt::SmartPtr<Ipopt::RegisteredOptions> RegOptions(){
    return roptions_;}
  
  /// Get the options (for getting theur values).
  virtual Ipopt::SmartPtr<const Ipopt::OptionsList> Options() const{
    return ConstPtr(options_);}
  
  /// Get the options (for getting and setting their values).
  virtual Ipopt::SmartPtr<Ipopt::OptionsList> Options(){
    return options_;}

   /// Get the CpuTime of the last optimization.
   virtual double CPUTime(){
   return (Ipopt::IsValid(cached_)) ? cached_->cpuTime_: 0.;}

   /// Get the iteration count of the last optimization.
   virtual int IterationCount(){
   return (Ipopt::IsValid(cached_)) ? cached_->istat[1]:0;}

  /// turn off all output from the solver 
  virtual void turnOffOutput()
  {if(Ipopt::IsValid(cached_)) cached_->iprint = 0;}
  /// turn on all output from the solver
  virtual void turnOnOutput()
  {if(Ipopt::IsValid(cached_)) cached_->iprint = 3;}

  /// Get the solver name
  virtual std::string & solverName(){
    return solverName_;}

   /// Register this solver options into passed roptions
   void RegisterOptions(){
   RegisterOptions(roptions_);}

   /// Register this solver options into passed roptions
   virtual void RegisterOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions);
private:
  /** @name Private function members. */
  /** @{ */
  /** Perform optimization using data structure in cache. */
  TNLPSolver::ReturnStatus callOptimizer();
  /** @} */
  /** Register all the options for filter. */
  

  /** Journalist */
  Ipopt::SmartPtr<Ipopt::Journalist> journalist_;

  /** Options */
  Ipopt::SmartPtr<Ipopt::OptionsList> options_;

  /** Registered Options */
  Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions_;
  


  /** Cached information for reoptimizing. */
  struct cachedInfo : public Ipopt::ReferencedObject {
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
    real * g_;
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
     g_(NULL),
     cstype(NULL),
     maxiter(1000),
     istat(NULL),
     rstat(NULL),
     tnlp_(NULL),
     hStruct_(NULL),
     permutationJac_(NULL),
     permutationHess_(NULL),
     cpuTime_(0)
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
     g_(NULL),
     cstype(NULL),
     maxiter(1000),
     istat(NULL),
     rstat(NULL),
     tnlp_(NULL),
     hStruct_(NULL),
     permutationJac_(NULL),
     permutationHess_(NULL),
     cpuTime_(0)   {
      initialize(tnlp, options);
    }

    /** Fill data structures for filter with info from tnlp. */
    void initialize(const Ipopt::SmartPtr<Ipopt::TNLP> &tnlp,
		    Ipopt::SmartPtr<Ipopt::OptionsList>& options);

    /** Optimize problem described by cache with filter.*/
    void optimize();

    /** Destructor. */
    ~cachedInfo(){
      delete [] x;
      delete []c;
      delete [] bounds;
      delete [] s;
      delete [] a;
      delete [] la;
      delete [] ws;
      delete [] lws;
      delete [] lam;
      delete [] g_;
      delete [] cstype;      
      delete [] istat;
      delete [] rstat;
      delete [] permutationJac_;
      delete [] permutationHess_;
      delete [] hStruct_;
      tnlp_ = NULL;
    }
  };

  /** Cached information on last problem optimized for reoptimization. */
  Ipopt::SmartPtr<cachedInfo> cached_;

  //name of solver (Filter)
  static std::string  solverName_;  
};

}// end namespace Bonmin
#endif

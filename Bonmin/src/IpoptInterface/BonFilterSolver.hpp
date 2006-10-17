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

namespace Bonmin{
class FilterSolver : public TNLPSolver{
public:
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

  virtual void enableWarmStart()
  {if(Ipopt::IsValid(cached_)) cached_->ifail= -1;}
  virtual void disableWarmStart()
  {if(Ipopt::IsValid(cached_)) cached_->ifail= 0;}
   //@}

   ///Get a pointer to RegisteredOptions (generally used to add new ones)
   virtual Ipopt::SmartPtr<Ipopt::RegisteredOptions> RegOptions();

   /// Get the options (for getting theur values).
   virtual Ipopt::SmartPtr<const Ipopt::OptionsList> Options() const;

   /// Get the options (for getting and setting their values).
   virtual Ipopt::SmartPtr<Ipopt::OptionsList> Options();

   /// Get the CpuTime of the last optimization.
   virtual double CPUTime();

   /// Get the iteration count of the last optimization.
   virtual int IterationCount(); 

  /// turn off all output from the solver 
  virtual void turnOffOutput()
  {if(Ipopt::IsValid(cached_)) cached_->iprint = 0;}
  /// turn on all output from the solver
  virtual void turnOnOutput()
  {if(Ipopt::IsValid(cached_)) cached_->iprint = 3;}
private:
  /** @name Private function members. */
  /** @{ */
  /** Perform optimization using data structure in cache. */
  TNLPSolver::ReturnStatus callOptimizer();

  void registerAllOptions();
  
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


    /** Constructor.*/
    cachedInfo():
     n(-1),
     m(-1),
     kmax(-1),
     maxa(-1),
     maxf(-1),
     mlp(-1),
     maxWk(-1),
     maxiWk(-1),
     iprint(-1),
     nout(-1),
     ifail(-100),
     rho(0),
      x(NULL),
      c(NULL),
     f(1e100),
     fmin(-1e100),
      bounds(NULL),
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
     hStruct_(NULL)
    {}

    cachedInfo(const Ipopt::SmartPtr<Ipopt::TNLP> &tnlp,
	       Ipopt::SmartPtr<Ipopt::OptionsList>& options)
    {
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
      delete [] a;
      delete [] la;
      delete [] ws;
      delete [] lws;
      delete [] lam;
      delete [] cstype;      
      delete [] istat;
      delete [] rstat;
      delete [] hStruct_;
      tnlp_ = NULL;
    }
  };
  Ipopt::SmartPtr<cachedInfo> cached_;
};

}// end namespace Bonmin
#endif

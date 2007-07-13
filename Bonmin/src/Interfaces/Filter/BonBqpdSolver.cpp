// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Andreas Waechter, International Business Machines Corporation
//                    based on BonFilterSolver.cpp
//
// Date : 07/09/2007

#include "BonminConfig.h"

#include "BonBqpdSolver.hpp"
#include "BonFilterWarmStart.hpp"

#include "CoinTime.hpp"

typedef Bonmin::BqpdSolver::fint fint;
typedef Bonmin::BqpdSolver::real real;

extern "C" {
  void F77_FUNC(bqpd,BQPD)(fint* n, fint* m, fint* k, fint* kmax,
			   real* a, fint* la, real* x, real* bl, real* bu,
			   real* f, real* fmin, real* g, real* r, real* w,
			   real* e, fint* ls, real* alp, fint* lp, fint* mlp,
			   fint* peq, real* ws, fint* lws, fint* m0de,
			   fint* ifail, fint* info, fint* iprint, fint* nout);

  //Access to filter common blocks
  extern struct {
    fint kk,ll,kkk,lll,mxws,mxlws;
  } F77_FUNC(wsc,WSC);

  extern struct {
    real eps,tol,emin;
  } F77_FUNC(epsc,EPSC);

  extern struct {
    real sgnf;
    fint nrep,npiv,nres;
  } F77_FUNC(repc,REPC);

  extern struct {
    fint nup,nfreq;
  } F77_FUNC(refactorc,REFACTORC);

  extern struct {
    real vstep;
  } F77_FUNC(vstepc,VSTEPC);

  extern struct {
    fint phl, phr, phc;
  } F77_FUNC(hessc,HESSC);

  extern struct {
    fint scale_mode, phe;
  } F77_FUNC(scalec,SCALEC);
}

namespace Bonmin{

  std::string BqpdSolver::solverName_ = "Bqpd QP";

  void BqpdSolver::
  registerOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions){
    roptions->SetRegisteringCategory("Bqpd options");
  }

  BqpdSolver::BqpdSolver(bool createEmpty /* = false */)
    :
    cached_(NULL)
  {
    options_ = new Ipopt::OptionsList();
    if (createEmpty) return;

    journalist_= new Ipopt::Journalist();
    roptions_ = new Ipopt::RegisteredOptions();
  
    try{
      Ipopt::SmartPtr<Ipopt::Journal> stdout_journal =
	journalist_->AddFileJournal("console", "stdout", Ipopt::J_ITERSUMMARY);
  
      registerOptions();

      options_->SetJournalist(journalist_);
      options_->SetRegisteredOptions(roptions_);
    }
    catch (Ipopt::IpoptException &E){
      E.ReportException(*journalist_);
      throw E;
    }
    catch(std::bad_alloc){
      journalist_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN, "\n Not enough memory .... EXIT\n");
      throw -1;
    }
    catch(...){
      Ipopt::IpoptException E("Uncaught exception in BqpdSolver::BqpdSolver()",
			      "BonBqpdSolver.cpp",-1);
      throw E;
    }
  }

  BqpdSolver::BqpdSolver(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions,
			 Ipopt::SmartPtr<Ipopt::OptionsList> options,
			 Ipopt::SmartPtr<Ipopt::Journalist> journalist)
    :
    journalist_(journalist),
    options_(options),
    roptions_(roptions),
    cached_(NULL)
  {
  }
                        
  Ipopt::SmartPtr <TNLPSolver>
  BqpdSolver::clone(){
    Ipopt::SmartPtr<BqpdSolver> retval = new BqpdSolver(true);
    *retval->options_ = *options_; // Copy the options
    retval->roptions_ = roptions_; // only copy pointers of registered options
    retval->journalist_ = journalist_; // and journalist
    return GetRawPtr(retval);
  }

  BqpdSolver::~BqpdSolver(){
  }

  bool
  BqpdSolver::Initialize(std::string optFile){
    std::ifstream is;
    if (optFile != "") {
      try {
        is.open(optFile.c_str());
      }
      catch(std::bad_alloc) {
        journalist_->Printf(Ipopt::J_SUMMARY, Ipopt::J_MAIN, "\nEXIT: Not enough memory.\n");
        return false;
      }
      catch(...) {
        Ipopt::IpoptException E("Unknown Exception caught in ipopt", "Unknown File", -1);
        E.ReportException(*journalist_);
        return false;
      }
    }
    bool retval = Initialize(is);
    if (is) {
      is.close();
    }
    return retval;
  }

  bool
  BqpdSolver::Initialize(std::istream &is){
    if(is.good()){
      options_->ReadFromStream(*journalist_, is);
    }
    return true;
  }

  /// Solves a problem expressed as a TQP
  TNLPSolver::ReturnStatus 
  BqpdSolver::OptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP>& tnlp)
  {
    BranchingTQP* tqp = dynamic_cast<BranchingTQP*>(GetRawPtr(tnlp));
    if (!tqp) {
      Ipopt::IpoptException E("BqpdSolver called with object other than a BranchingTQP",
			      "BonBqpdSolver.cpp",-1);
      throw E;
    }
    if (IsNull(cached_) || !cached_->use_warm_start_in_cache_) {
      cached_ = new cachedInfo(tqp, options_);
    }
    return callOptimizer();
  }

/// Solves a problem expresses as a TNLP 
  TNLPSolver::ReturnStatus 
  BqpdSolver::ReOptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> & tnlp)
  {
#ifndef NDEBUG
    BranchingTQP* tqp = dynamic_cast<BranchingTQP*>(GetRawPtr(tnlp));
    assert(tqp == GetRawPtr(cached_->tqp_));
#endif
    int n = cached_->n;
    int m = cached_->m;
    tnlp->get_bounds_info(n, cached_->bl, cached_->bu, 
			  m, cached_->bl+n, cached_->bu+n);

    return callOptimizer();
  }



void
BqpdSolver::cachedInfo::initialize(const Ipopt::SmartPtr<BranchingTQP> & tqp,
				   Ipopt::SmartPtr<Ipopt::OptionsList>& options)
{
  // In case BQPD's BLOCK DATA doesn't work, we initialize the COMMON
  // BLOCKs explicitly here
  F77_FUNC(epsc,EPSC).eps = 1111e-19;
  F77_FUNC(epsc,EPSC).tol = 1e-12;
  F77_FUNC(epsc,EPSC).emin = 1.;
  F77_FUNC(repc,REPC).sgnf = 1e-4;
  F77_FUNC(repc,REPC).nrep = 2;
  F77_FUNC(repc,REPC).npiv = 3;
  F77_FUNC(repc,REPC).nres = 2;
  F77_FUNC(refactorc,REFACTORC).nfreq = 500;

  Ipopt::TNLP::IndexStyleEnum index_style;
  Index nv, nc, nnz_jac_g, nnz_hess;
  tqp->get_nlp_info(nv, nc, nnz_jac_g, nnz_hess, index_style);
  n = nv;
  m = nc;

  Index kmax_ipt;
  options->GetIntegerValue("kmax", kmax_ipt, "bqpd.");
  if (kmax_ipt == -1) {
    kmax = n;
  }
  else {
    kmax = kmax_ipt;
    kmax = min(kmax,n);
  }
  Index mlp_ipt;
  options->GetIntegerValue("mlp", mlp_ipt,"bqpd.");
  mlp = mlp_ipt;

  Index mxws_ipt;
  options->GetIntegerValue("mxws", mxws_ipt, "bqpd.");
  mxws = mxws_ipt;

  Ipopt::Index mxlws_ipt;
  options->GetIntegerValue("mxlws",  mxlws_ipt, "bqpd.");
  mxlws = mxlws_ipt;

  use_warm_start_in_cache_ = false;

  // Get space for arrays
  x = new real[n];
  g = new real[n];
  bl = new real[n+m];
  bu = new real[n+m];
  g = new real[n];
  r = new real[n+m];
  w = new real[n+m];
  e = new real[n+m];
  ls = new fint[n+m];
  alp = new real[mlp];
  lp = new fint[mlp];
  ws = new real[mxws];
  lws = new fint[mxlws];

  // Getting the bounds
  tqp->get_bounds_info(n, bl, bu, m, bl+n, bu+n);
  const double infty = 1e30;
  for(int i = 0 ; i < n+m ; i++){
    if(bl[i] < -infty) bl[i] = - infty;
  }
  for(int i = 0 ; i < n+m ; i++){
    if(bu[i] > infty) bu[i] = infty;
  }

  // Set up sparse matrix with objective gradient and constraint Jacobian

  const Number* obj_grad = tqp->ObjGrad();
  int amax = nnz_jac_g;
  for(int i = 0; i<n; i++) {
    if (obj_grad[i]!=0.) {
      amax++;
    }
  }
  int lamax = amax+m+2;
  a = new real[amax];
  la = new fint [1+lamax];

  // Objective function gradient
  int nnz_grad = 0;
  for(int i = 1; i < n ; i++) {
    if (obj_grad[i]!=0.) {
      a[nnz_grad] = obj_grad[i];
      la[++nnz_grad] = i;
    }
  }

  // Constraint Jacobian
  const Number* JacVals = tqp->ConstrJacVals();
  const Index* RowJac = tqp->ConstrJacIRow();
  const Index* ColJac = tqp->ConstrJacJCol();

  int* permutationJac = new int [nnz_jac_g];
  FilterSolver::TMat2RowPMat(n, m, nnz_jac_g,  RowJac, ColJac, permutationJac,
			     la, nnz_grad);
  for(int i=0; i<nnz_jac_g; i++) {
    const int& indice = permutationJac[i];
    a[nnz_grad+i] = JacVals[indice];
  }
  delete [] permutationJac;
  
  // Now setup Hessian
  const Number* HessVals = tqp->ObjHessVals();
  const Index* RowHess = tqp->ObjHessIRow();
  const Index* ColHess = tqp->ObjHessJCol();

  kk = nnz_hess;
  ll = nnz_hess + n + 2;
  int* permutationHess = new int[nnz_hess];

  FilterSolver::TMat2ColPMat(n, m, nnz_hess, RowHess, ColHess,
			     permutationHess, lws, 0);
  for(int i=0; i<nnz_hess; i++) {
    ws[i] = HessVals[permutationHess[i]];
  }
  delete [] permutationHess;

  Index bufy;
  options->GetIntegerValue("iprint",bufy, "bqpd.");
  iprint = bufy;
  nout = 6;

  tqp_ = tqp;
}

/// Solves a problem expresses as a TNLP 
TNLPSolver::ReturnStatus 
BqpdSolver::callOptimizer()
{
  cached_->optimize();

  TNLPSolver::ReturnStatus optimizationStatus = TNLPSolver::exception;
  Ipopt::SolverReturn status = Ipopt::INTERNAL_ERROR;
  fint ifail = cached_->ifail;
  switch(ifail){
  case 0:
    optimizationStatus = TNLPSolver::solvedOptimal;
    status = Ipopt::SUCCESS;
    break;
  case 1:
    optimizationStatus = TNLPSolver::unbounded;
    status = Ipopt::DIVERGING_ITERATES;
  case 2:
  case 3:
    optimizationStatus = TNLPSolver::provenInfeasible;
    status = Ipopt::LOCAL_INFEASIBILITY;
    break;
  }

  Index dummy_len = Ipopt::Max(cached_->n,cached_->m);
  Number* dummy = new Number[dummy_len];
  for(int i=0; i<dummy_len; i++) {
    dummy[i] = 0.;
  }

  cached_->tqp_->finalize_solution(status, cached_->n, 
				   cached_->x, dummy, dummy, 
				   cached_->m, dummy, dummy, 
				   cached_->f, NULL, NULL);
  return optimizationStatus;
}

/** Optimize problem described by cache with filter.*/
void 
BqpdSolver::cachedInfo::optimize()
{
  printf("blabla");
  if (use_warm_start_in_cache_) {
    m0de = 6;
    use_warm_start_in_cache_ = false;
  }
  else {
    m0de = 0;
    tqp_->get_starting_point(n, 1, x, 0, NULL, NULL, m, 0, NULL);
    ifail = 0;
  }

  // Set up some common block stuff
  F77_FUNC(scalec,SCALEC).scale_mode = 0;  // No scaling
  F77_FUNC(scalec,SCALEC).phe = 0;  // No scaling
  F77_FUNC(hessc,HESSC).phl = 1; // This is to tell gdotx to do the right thing
  F77_FUNC(wsc,WSC).kk = kk;
  F77_FUNC(wsc,WSC).ll = ll;
  F77_FUNC(wsc,WSC).mxws = mxws;
  F77_FUNC(wsc,WSC).mxlws = mxlws;

  printf("mode = %d vstep = %e tol = %e\n", m0de, F77_FUNC(vstepc,VSTEPC).vstep,F77_FUNC(epsc,EPSC).tol);

  cpuTime_ = - CoinCpuTime();
  real fmin = -1e100;
  F77_FUNC(bqpd,BQPD)(&n, &m, &k, &kmax, a, la, x, bl, bu, &f, &fmin,
		      g, r, w, e, ls, alp, lp, &mlp, &peq, ws, lws,
		      &m0de, &ifail, info, &iprint, &nout);
  cpuTime_ += CoinCpuTime();
}

std::string
BqpdSolver::UnsolvedBqpdError::errorNames_[1] =
{"Internal error in Filter SQP."};

std::string
BqpdSolver::UnsolvedBqpdError::solverName_ = 
          "filterSqp";

const std::string& 
BqpdSolver::UnsolvedBqpdError::errorName() const
{ return errorNames_[0];}

const std::string& 
BqpdSolver::UnsolvedBqpdError::solverName() const
{ return solverName_;}

bool
BqpdSolver::setWarmStart(const CoinWarmStart * warm, 
			   Ipopt::SmartPtr<TMINLP2TNLP> tnlp)
{
#if 0
  if (IsNull(cached_)) {
    cached_ = new cachedInfo(GetRawPtr(tnlp), options_);
  }

  const FilterWarmStart * warmF = dynamic_cast<const FilterWarmStart *> (warm);
  //CoinCopyN(warmF->xArray(), warmF->xSize(), cached_->x);
  const fint xsize = warmF->xSize();
  real* x = cached_->x;
  const real* xarray = warmF->xArray();
  for (int i = 0; i<xsize; i++) {
    x[i] = xarray[i];
  }
  CoinCopyN(warmF->lamArray(), warmF->lamSize(), cached_->lam);
  CoinCopyN(warmF->lwsArray(), warmF->lwsSize(), cached_->lws);
  for(int i = 0 ; i < 14 ; i ++) {
    cached_->istat[i] = warmF->istat()[i];
  }
  cached_->use_warm_start_in_cache_ = true;
#endif
  printf("BqpdSolver::setWarmStart called!\n");
  return true;
}

CoinWarmStart *
BqpdSolver::getWarmStart(Ipopt::SmartPtr<TMINLP2TNLP> tnlp) const
{
#if 0
  return new FilterWarmStart(cached_->n, cached_->x,
			     cached_->n+cached_->m, cached_->lam,
			     cached_->maxiWk, cached_->lws, cached_->istat);
#endif
  printf("BqpdSolver::getWarmStart called!\n");
  return NULL;
}

CoinWarmStart * 
BqpdSolver::getEmptyWarmStart() const {
#if 0
  return new FilterWarmStart;
#endif
  printf("BqpdSolver::getEmptyWarmStart called \n");
  return NULL;
}

}//end namespace Bonmin

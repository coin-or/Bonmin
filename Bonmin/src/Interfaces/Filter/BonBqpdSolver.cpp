// (C) Copyright International Business Machines Corporation 2007, 2008
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Andreas Waechter, International Business Machines Corporation
//                    based on BonFilterSolver.cpp
//
// Date : 07/09/2007

#include "BonminConfig.h"

#include "BonBqpdSolver.hpp"
#include "BonBqpdWarmStart.hpp"

#include "CoinTime.hpp"
#include <algorithm>

#define InitializeAll

typedef Bonmin::BqpdSolver::fint fint;
typedef Bonmin::BqpdSolver::real real;
int Bonmin::BqpdSolver::reinit_freq_ = 0;
int Bonmin::BqpdSolver::m0de_ = 6;
extern "C"
{
  void F77_FUNC(bqpd,BQPD)(fint* n, fint* m, fint* k, fint* kmax,
      real* a, fint* la, real* x, real* bl, real* bu,
      real* f, real* fmin, real* g, real* r, real* w,
      real* e, fint* ls, real* alp, fint* lp, fint* mlp,
      fint* peq, real* ws, fint* lws, fint* m0de,
      fint* ifail, fint* info, fint* iprint, fint* nout);

  //Access to filter common blocks
  extern struct {
      fint kk,ll,kkk,lll,mxws,mxlws;
    }
  F77_FUNC(wsc,WSC);

  extern struct {
      real eps,tol,emin;
    }
  F77_FUNC(epsc,EPSC);

  extern struct {
      real sgnf;
      fint nrep,npiv,nres;
    }
  F77_FUNC(repc,REPC);

  extern struct {
      fint nup,nfreq;
    }
  F77_FUNC(refactorc,REFACTORC);

  extern struct {
      real vstep;
    }
  F77_FUNC(vstepc,VSTEPC);

  extern struct {
      fint phl, phr, phc;
    }
  F77_FUNC(hessc,HESSC);

  extern struct {
      fint scale_mode, phe;
    }
  F77_FUNC(scalec,SCALEC);

  extern struct {
      fint irh1,na,na1,nb,nb1,ka1,kb1,kc1,irg1,lu1,lv,lv1,ll1;
    }
  F77_FUNC(bqpdc,BQPDC);

  extern struct {
      real alpha;
    }
  F77_FUNC(alphac,ALPHAC);

  extern struct {
      fint ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof,lc;
      fint lc1,li,li1,lm,lm1,lp_,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1;
    }
  F77_FUNC(sparsec,SPARSEC);

  extern struct {
      fint m1,m2,mp,mq,lastr,irow;
    }
  F77_FUNC(factorc,FACTORC);

  extern struct {
      fint mxm1;
    }
  F77_FUNC(mxm1c,MXM1C);

  extern struct {
      real c;
    }
  F77_FUNC(minorc,MINORS);
}

namespace Bonmin
{

  std::string BqpdSolver::solverName_ = "Bqpd QP";

  void BqpdSolver::
  registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("Bqpd options",RegisteredOptions::BqpdCategory);
    roptions->AddLowerBoundedNumberOption("qp_fillin_factor", "Factor for estimating fill-in for factorization method in Bqpd", 0., true, 4.);
    roptions->AddBoundedIntegerOption("hot_start_m0de", "Choice of the hot start option", 0, 6, 6);
    roptions->AddLowerBoundedIntegerOption("reinit_freq", "Number of pivots before recopy hot start", 0, 0);
  }

  BqpdSolver::BqpdSolver(bool createEmpty /* = false */)
      :
      TNLPSolver(),
      cached_(NULL)
  {
    if (createEmpty) return;
  }

  BqpdSolver::BqpdSolver(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
      Ipopt::SmartPtr<Ipopt::OptionsList> options,
      Ipopt::SmartPtr<Ipopt::Journalist> journalist,
      const std::string & prefix)
      :
      TNLPSolver(roptions, options, journalist, prefix),
      cached_(NULL)
  {
    options->GetNumericValue("qp_fillin_factor", fillin_factor_, "bqpd.");
    options->GetIntegerValue("kmax", kmax_ipt_, "bqpd.");
    options->GetIntegerValue("mlp", mlp_ipt_,"bqpd.");
    options->GetIntegerValue("hot_start_m0de", m0de_,"bqpd.");
    options->GetIntegerValue("reinit_freq", reinit_freq_,"bqpd.");
     if(!options_->GetIntegerValue("print_level",default_log_level_,""))
      default_log_level_ = 1;

  }

  Ipopt::SmartPtr <TNLPSolver>
  BqpdSolver::clone()
  {
#ifdef TIME_BQPD
    times().create -= CoinCpuTime();
#endif
    Ipopt::SmartPtr<BqpdSolver> retval = new BqpdSolver(true);
    *retval->options_ = *options_; // Copy the options
    retval->roptions_ = roptions_; // only copy pointers of registered options
    retval->journalist_ = journalist_; // and journalist
    retval->prefix_ = prefix_;
    retval->fillin_factor_ = fillin_factor_;
    retval->kmax_ipt_ = kmax_ipt_;
    retval->mlp_ipt_ = mlp_ipt_;
    retval->default_log_level_ = default_log_level_;
    retval->reinit_freq_ = reinit_freq_;
    retval->m0de_ = m0de_;
#ifdef TIME_BQPD
    times().create += CoinCpuTime();
#endif
    return GetRawPtr(retval);
  }

  BqpdSolver::~BqpdSolver()
  {}

  bool
  BqpdSolver::Initialize(std::string optFile)
  {
    std::ifstream is;
    if (optFile != "") {
      try {
        is.open(optFile.c_str());
      }
      catch (std::bad_alloc) {
        journalist_->Printf(Ipopt::J_SUMMARY, Ipopt::J_MAIN, "\nEXIT: Not enough memory.\n");
        return false;
      }
#ifndef NO_CATCH_ALL
      catch (...) {
        Ipopt::IpoptException E("Unknown Exception caught in ipopt", "Unknown File", -1);
        E.ReportException(*journalist_);
        return false;
      }
#endif
    }
    bool retval = Initialize(is);
    if (is) {
      is.close();
    }
    if(!options_->GetIntegerValue("print_level",default_log_level_,""))
      default_log_level_ = 1;
    return retval;
  }

  bool
  BqpdSolver::Initialize(std::istream &is)
  {
    if (is.good()) {
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
      cached_ = new cachedInfo(tqp, options_, kmax_ipt_, mlp_ipt_,
			       &fillin_factor_);
    }
#ifdef TIME_BQPD
    times().solve -= CoinCpuTime();
#endif
    TNLPSolver::ReturnStatus r_val = callOptimizer();
#ifdef TIME_BQPD
    times().solve += CoinCpuTime();
#endif
    return r_val;
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
    // Make sure bounds are not infinity
    for (int i=0; i<n+m; i++) {
      cached_->bl[i] = std::max(cached_->bl[i], -1e50);
      cached_->bu[i] = std::min(cached_->bu[i], 1e50);
    }

#if 1
    cached_->use_warm_start_in_cache_ = true;  // Trying...
    cached_->m0de = m0de_;
#else
    cached_->re_initialize();
    cached_->use_warm_start_in_cache_ = true;  // Trying...
#endif

#ifdef TIME_BQPD
    times().resolve -= CoinCpuTime();
#endif
    TNLPSolver::ReturnStatus r_val = callOptimizer();
#ifdef TIME_BQPD
    times().resolve += CoinCpuTime();
#endif
    return r_val;
  }

  void
  BqpdSolver::cachedInfo::initialize(const Ipopt::SmartPtr<BranchingTQP> & tqp,
				     Ipopt::SmartPtr<Ipopt::OptionsList>& options,
				     int kmax_ipt,
				     int mlp_ipt,
				     double* fillin_factor)
  {
    // Maybe a dirty trick?  We want to change fillin_factor in calling object
    fillin_factor_ = fillin_factor;

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
    Ipopt::Index nv, nc, nnz_jac_g, nnz_hess;
    tqp->get_nlp_info(nv, nc, nnz_jac_g, nnz_hess, index_style);
    n = nv;
    m = nc;

    if (kmax_ipt == -1) {
      kmax = n;
    }
    else {
      kmax = kmax_ipt;
      kmax = std::min(kmax,n);
    }
    mlp = mlp_ipt;

    use_warm_start_in_cache_ = false;

    // Get space for arrays
    x = new real[n];
    bl = new real[n+m];
    bu = new real[n+m];
    g = new real[n];
    r = new real[n+m];
    w = new real[n+m];
    e = new real[n+m];
    ls = new fint[n+m];
    alp = new real[mlp];
    lp = new fint[mlp];

    // Getting the bounds
    tqp->get_bounds_info(n, bl, bu, m, bl+n, bu+n);

    // Make sure bounds are not infinity
    for (int i=0; i<n+m; i++) {
      bl[i] = std::max(bl[i], -1e50);
      bu[i] = std::min(bu[i], 1e50);
    }

    // Set up sparse matrix with objective gradient and constraint Jacobian

    const Ipopt::Number* obj_grad = tqp->ObjGrad();
    amax_ = nnz_jac_g;
    for (int i = 0; i<n; i++) {
      if (obj_grad[i]!=0.) {
        amax_++;
      }
    }
    int lamax = amax_+m+2;
    a = new real[amax_];
    la = new fint [1+lamax];

    // Objective function gradient
    int nnz_grad = 0;
    for (int i = 0; i < n ; i++) {
      if (obj_grad[i]!=0.) {
        a[nnz_grad] = obj_grad[i];
        la[++nnz_grad] = i+1;
      }
    }
    la[amax_+1] = 1;

    // Constraint Jacobian
    const Ipopt::Number* JacVals = tqp->ConstrJacVals();
    const Ipopt::Index* RowJac = tqp->ConstrJacIRow();
    const Ipopt::Index* ColJac = tqp->ConstrJacJCol();

    int* permutationJac = new int [nnz_jac_g];
    FilterSolver::TMat2RowPMat(false, n, m, nnz_jac_g,  RowJac, ColJac, permutationJac,
        la, nnz_grad, 1, Ipopt::TNLP::C_STYLE);
    for (int i=0; i<nnz_jac_g; i++) {
      const int& indice = permutationJac[i];
      a[nnz_grad+i] = JacVals[indice];
    }
    delete [] permutationJac;
#if 0
//deleteme
    printf("nnz_grad = %d nnz_jac = %d\n", nnz_grad, nnz_jac_g);
    for (int i=0; i<1+lamax; i++) printf("la[%2d] = %d\n", i,la[i]);
    for (int i=0; i<amax_; i++) printf("a[%3d] = %e\n",i,a[i]);
#endif

    // Setup workspaces
    mxws = nnz_hess + ((kmax*(kmax+9))/2 + 2*n + m) + (5*n + (int)(*fillin_factor_*(double)amax_));
    mxlws = (nnz_hess + n + 2) + kmax + (9*n + m);

    ws = new real[mxws];
    lws = new fint[mxlws];

#ifdef InitializeAll
    for (int i=0; i<mxws; i++) {
      ws[i] = 42.;
    }
    for (int i=0; i<mxlws; i++) {
      lws[i] = 55;
    }
#endif

    // Now setup Hessian
    const Ipopt::Number* HessVals = tqp->ObjHessVals();
    const Ipopt::Index* RowHess = tqp->ObjHessIRow();
    const Ipopt::Index* ColHess = tqp->ObjHessJCol();

    kk = nnz_hess;
    ll = nnz_hess + n + 2;
    int* permutationHess = new int[nnz_hess];

    FilterSolver::TMat2RowPMat(true, n, n, nnz_hess, RowHess, ColHess,
        permutationHess, lws, 0, 0, Ipopt::TNLP::C_STYLE);
    for (int i=0; i<nnz_hess; i++) {
      ws[i] = HessVals[permutationHess[i]];
    }
    delete [] permutationHess;
#if 0
//deleteme
    printf("nnz_hess = %d\n", nnz_hess);
    for (int i=0; i<ll; i++) printf("hess lws[%2d] = %d\n", i,lws[i]);
    for (int i=0; i<kk; i++) printf("hess ws[%3d] = %e\n",i,ws[i]);
#endif

    Ipopt::Index bufy;
    options->GetIntegerValue("iprint",bufy, "bqpd.");
    iprint = bufy;
    nout = 6;

    tqp_ = tqp;
  }

  void
  BqpdSolver::cachedInfo::re_initialize()
  {
    use_warm_start_in_cache_ = false;

    // Setup workspaces
    for (int i=0; i<mxws; i++) {
      ws[i] = 42.;
    }
    for (int i=0; i<mxlws; i++) {
      lws[i] = 55;
    }
  }

/// Solves a problem expresses as a TNLP
  TNLPSolver::ReturnStatus BqpdSolver::callOptimizer() {
#ifdef TIME_BQPD
    times().numsolve++;
#endif
    cached_->optimize();

    TNLPSolver::ReturnStatus optimizationStatus = TNLPSolver::exception;
    Ipopt::SolverReturn status = Ipopt::INTERNAL_ERROR;
    fint ifail = cached_->ifail;
    switch (ifail) {
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

    Ipopt::Index dummy_len = std::max(cached_->n,cached_->m);
    Ipopt::Number* dummy = new Ipopt::Number[dummy_len];
    for (int i=0; i<dummy_len; i++) {
      dummy[i] = 0.;
    }

    cached_->tqp_->finalize_solution(status, cached_->n,
        cached_->x, dummy, dummy,
        cached_->m, dummy, dummy,
        cached_->f, NULL, NULL);
    delete [] dummy;
    return optimizationStatus;
  }

  BqpdSolver::cachedInfo::~cachedInfo()
  {
    if (haveHotStart_) {
      unmarkHotStart();
    }
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

  bool
  BqpdSolver::cachedInfo::markHotStart()
  {
     next_reinit_ = BqpdSolver::reinit_freq_; 
#ifdef DISABLE_COPYING
     return 1;
#endif
#ifdef TIME_BQPD
    times_.warm_start -= CoinCpuTime();
#endif
    haveHotStart_ = true;
    irh1 = F77_FUNC(bqpdc,BQPDC).irh1;
    na = F77_FUNC(bqpdc,BQPDC).na;
    na1 = F77_FUNC(bqpdc,BQPDC).na1;
    nb = F77_FUNC(bqpdc,BQPDC).nb;
    nb1 = F77_FUNC(bqpdc,BQPDC).nb1;
    ka1 = F77_FUNC(bqpdc,BQPDC).ka1;
    kb1 = F77_FUNC(bqpdc,BQPDC).kb1;
    kc1 = F77_FUNC(bqpdc,BQPDC).kc1;
    irg1 = F77_FUNC(bqpdc,BQPDC).irg1;
    lu1 = F77_FUNC(bqpdc,BQPDC).lu1;
    lv = F77_FUNC(bqpdc,BQPDC).lv;
    lv1 = F77_FUNC(bqpdc,BQPDC).lv1;
    ll1 = F77_FUNC(bqpdc,BQPDC).ll1;
    eps = F77_FUNC(epsc,EPSC).eps;
    tol = F77_FUNC(epsc,EPSC).tol;
    emin = F77_FUNC(epsc,EPSC).emin;
    vstep = F77_FUNC(vstepc,VSTEPC).vstep;
    sgnf = F77_FUNC(repc,REPC).sgnf;
    nrep = F77_FUNC(repc,REPC).nrep;
    npiv = F77_FUNC(repc,REPC).npiv;
    nres = F77_FUNC(repc,REPC).nres;
    nup = F77_FUNC(refactorc,REFACTORC).nup;
    nfreq = F77_FUNC(refactorc,REFACTORC).nfreq;
    alpha = F77_FUNC(alphac,ALPHAC).alpha;

    ns = F77_FUNC(sparsec,SPARSEC).ns;
    ns1 = F77_FUNC(sparsec,SPARSEC).ns1;
    nt = F77_FUNC(sparsec,SPARSEC).nt;
    nt1 = F77_FUNC(sparsec,SPARSEC).nt1;
    nu = F77_FUNC(sparsec,SPARSEC).nu;
    nu1 = F77_FUNC(sparsec,SPARSEC).nu1;
    nx = F77_FUNC(sparsec,SPARSEC).nx;
    nx1 = F77_FUNC(sparsec,SPARSEC).nx1;
    np = F77_FUNC(sparsec,SPARSEC).np;
    np1 = F77_FUNC(sparsec,SPARSEC).np1;
    nprof = F77_FUNC(sparsec,SPARSEC).nprof;
    lc = F77_FUNC(sparsec,SPARSEC).lc;
    lc1 = F77_FUNC(sparsec,SPARSEC).lc1;
    li = F77_FUNC(sparsec,SPARSEC).li;
    li1 = F77_FUNC(sparsec,SPARSEC).li1;
    lm = F77_FUNC(sparsec,SPARSEC).lm;
    lm1 = F77_FUNC(sparsec,SPARSEC).lm1;
    lp_ = F77_FUNC(sparsec,SPARSEC).lp_;
    lp1 = F77_FUNC(sparsec,SPARSEC).lp1;
    lq = F77_FUNC(sparsec,SPARSEC).lq;
    lq1 = F77_FUNC(sparsec,SPARSEC).lq1;
    lr = F77_FUNC(sparsec,SPARSEC).lr;
    lr1 = F77_FUNC(sparsec,SPARSEC).lr1;
    ls_ = F77_FUNC(sparsec,SPARSEC).ls_;
    ls1 = F77_FUNC(sparsec,SPARSEC).ls1;
    lt = F77_FUNC(sparsec,SPARSEC).lt;
    lt1 = F77_FUNC(sparsec,SPARSEC).lt1;

    m1 = F77_FUNC(factorc,FACTORC).m1;
    m2 = F77_FUNC(factorc,FACTORC).m2;
    mp = F77_FUNC(factorc,FACTORC).mp;
    mq = F77_FUNC(factorc,FACTORC).mq;
    lastr = F77_FUNC(factorc,FACTORC).lastr;
    irow = F77_FUNC(factorc,FACTORC).irow;

    mxm1 = F77_FUNC(mxm1c,MXM1c).mxm1;
    c = F77_FUNC(minorc,MINORC).c;

    kkHot = F77_FUNC(wsc,WSC).kk;
    kkkHot = F77_FUNC(wsc,WSC).kkk;
    llHot = F77_FUNC(wsc,WSC).ll;
    lllHot = F77_FUNC(wsc,WSC).lll;

    kHot = k;
    xHot = CoinCopyOfArray(x,n);
    fHot = f;
    gHot = CoinCopyOfArray(g,n);
    rHot = CoinCopyOfArray(r,n+m);
    wHot = CoinCopyOfArray(w,n+m);
    eHot = CoinCopyOfArray(e,n+m);
    lsHot = CoinCopyOfArray(ls,n+m);
    alpHot = CoinCopyOfArray(alp,mlp);
    lpHot = CoinCopyOfArray(lp,mlp);
    peqHot = peq;
    wsHot = CoinCopyOfArray(ws,mxws);
    lwsHot = CoinCopyOfArray(lws,mxlws);
    infoHot[0] = info[0];
#ifdef TIME_BQPD
    times_.warm_start += CoinCpuTime();
#endif

    return true;
  }

  void
  BqpdSolver::cachedInfo::copyFromHotStart()
  {
#ifdef DISABLE_COPYING
    return;
#endif
#ifdef TIME_BQPD
    times_.warm_start -= CoinCpuTime();
#endif
    F77_FUNC(bqpdc,BQPDC).irh1 = irh1;
    F77_FUNC(bqpdc,BQPDC).na = na;
    F77_FUNC(bqpdc,BQPDC).na1 = na1;
    F77_FUNC(bqpdc,BQPDC).nb = nb;
    F77_FUNC(bqpdc,BQPDC).nb1 = nb1;
    F77_FUNC(bqpdc,BQPDC).ka1 = ka1;
    F77_FUNC(bqpdc,BQPDC).kb1 = kb1;
    F77_FUNC(bqpdc,BQPDC).kc1 = kc1;
    F77_FUNC(bqpdc,BQPDC).irg1 = irg1;
    F77_FUNC(bqpdc,BQPDC).lu1 = lu1;
    F77_FUNC(bqpdc,BQPDC).lv = lv;
    F77_FUNC(bqpdc,BQPDC).lv1 = lv1;
    F77_FUNC(bqpdc,BQPDC).ll1 = ll1;
    F77_FUNC(epsc,EPSC).eps = eps;
    F77_FUNC(epsc,EPSC).tol = tol;
    F77_FUNC(epsc,EPSC).emin = emin;
    F77_FUNC(vstepc,VSTEPC).vstep = vstep;
    F77_FUNC(repc,REPC).sgnf = sgnf;
    F77_FUNC(repc,REPC).nrep = nrep;
    F77_FUNC(repc,REPC).npiv = npiv;
    F77_FUNC(repc,REPC).nres = nres;
    F77_FUNC(refactorc,REFACTORC).nup = nup;
    F77_FUNC(refactorc,REFACTORC).nfreq = nfreq;
    F77_FUNC(alphac,ALPHAC).alpha = alpha;

    F77_FUNC(sparsec,SPARSEC).ns = ns;
    F77_FUNC(sparsec,SPARSEC).ns1 = ns1;
    F77_FUNC(sparsec,SPARSEC).nt = nt;
    F77_FUNC(sparsec,SPARSEC).nt1 = nt1;
    F77_FUNC(sparsec,SPARSEC).nu = nu;
    F77_FUNC(sparsec,SPARSEC).nu1 = nu1;
    F77_FUNC(sparsec,SPARSEC).nx = nx;
    F77_FUNC(sparsec,SPARSEC).nx1 = nx1;
    F77_FUNC(sparsec,SPARSEC).np = np;
    F77_FUNC(sparsec,SPARSEC).np1 = np1;
    F77_FUNC(sparsec,SPARSEC).nprof = nprof;
    F77_FUNC(sparsec,SPARSEC).lc = lc;
    F77_FUNC(sparsec,SPARSEC).lc1 = lc1;
    F77_FUNC(sparsec,SPARSEC).li = li;
    F77_FUNC(sparsec,SPARSEC).li1 = li1;
    F77_FUNC(sparsec,SPARSEC).lm = lm;
    F77_FUNC(sparsec,SPARSEC).lm1 = lm1;
    F77_FUNC(sparsec,SPARSEC).lp_ = lp_;
    F77_FUNC(sparsec,SPARSEC).lp1 = lp1;
    F77_FUNC(sparsec,SPARSEC).lq = lq;
    F77_FUNC(sparsec,SPARSEC).lq1 = lq1;
    F77_FUNC(sparsec,SPARSEC).lr = lr;
    F77_FUNC(sparsec,SPARSEC).lr1 = lr1;
    F77_FUNC(sparsec,SPARSEC).ls_ = ls_;
    F77_FUNC(sparsec,SPARSEC).ls1 = ls1;
    F77_FUNC(sparsec,SPARSEC).lt = lt;
    F77_FUNC(sparsec,SPARSEC).lt1 = lt1;

    F77_FUNC(factorc,FACTORC).m1 = m1;
    F77_FUNC(factorc,FACTORC).m2 = m2;
    F77_FUNC(factorc,FACTORC).mp = mp;
    F77_FUNC(factorc,FACTORC).mq = mq;
    F77_FUNC(factorc,FACTORC).lastr = lastr;
    F77_FUNC(factorc,FACTORC).irow = irow;

    F77_FUNC(mxm1c,MXM1c).mxm1 = mxm1;
    F77_FUNC(minorc,MINORC).c = c;

    F77_FUNC(wsc,WSC).kk = kkHot;
    F77_FUNC(wsc,WSC).kkk = kkkHot;
    F77_FUNC(wsc,WSC).ll = llHot;
    F77_FUNC(wsc,WSC).lll = lllHot;

    k = kHot;
    CoinCopyN(xHot,n,x);
    f = fHot;
    CoinCopyN(gHot,n,g);
    CoinCopyN(rHot,n+m,r);
    CoinCopyN(wHot,n+m,w);
    CoinCopyN(eHot,n+m,e);
    CoinCopyN(lsHot,n+m,ls);
    CoinCopyN(alpHot,mlp,alp);
    CoinCopyN(lpHot,mlp,lp);
    CoinCopyN(wsHot,mxws,ws);
    CoinCopyN(lwsHot,mxlws,lws);
    info[0] = infoHot[0];
    peq = peqHot;
#ifdef TIME_BQPD
    times_.warm_start += CoinCpuTime();
#endif

  }

  void
  BqpdSolver::cachedInfo::unmarkHotStart()
  {
    delete [] xHot;
    delete [] gHot;
    delete [] rHot;
    delete [] wHot;
    delete [] eHot;
    delete [] lsHot;
    delete [] alpHot;
    delete [] lpHot;
    delete [] wsHot;
    delete [] lwsHot;
    haveHotStart_ = false;
  }

  /** Optimize problem described by cache with Bqpd.*/
  void
  BqpdSolver::cachedInfo::optimize()
  {
    // Set up some common block stuff
    F77_FUNC(scalec,SCALEC).scale_mode = 0;  // No scaling
    F77_FUNC(scalec,SCALEC).phe = 0;  // No scaling
    F77_FUNC(hessc,HESSC).phl = 1; // This is to tell gdotx to do the right thing
    F77_FUNC(wsc,WSC).kk = kk;
    F77_FUNC(wsc,WSC).ll = ll;
    F77_FUNC(wsc,WSC).mxws = mxws;
    F77_FUNC(wsc,WSC).mxlws = mxlws;

    if (use_warm_start_in_cache_ && !bad_warm_start_info_) {
      ifail = 0;
      use_warm_start_in_cache_ = false;
      if (haveHotStart_ && pivots_ > next_reinit_) {
        //printf("Reinitialize hot start\n");
        copyFromHotStart();
        while (BqpdSolver::reinit_freq_ > 0&& next_reinit_ < pivots_)
          next_reinit_ += BqpdSolver::reinit_freq_;
      }
    }
    else {
      m0de = 0;
      tqp_->get_starting_point(n, 1, x, 0, NULL, NULL, m, 0, NULL);
      ifail = 0;
      bad_warm_start_info_ = false;
    }

#if 0
    printf("========= 222222222222 =============\n");
    printf("kk = %d ll = %d mxws = %d mxlws = %d\n", kk, ll, mxws, mxlws);
    for (int i=0; i<n; i++) {
      printf("xL[%3d] = %15.8e  xU[%3d] = %15.8e\n", i, bl[i], i, bu[i]);
    }
    for (int i=0; i<m; i++) {
      printf("gL[%3d] = %15.8e  gU[%3d] = %15.8e\n", i, bl[n+i], i, bu[n+i]);
    }
#endif
    cpuTime_ = - CoinCpuTime();
    real fmin = -1e100;
#if 0
    for (int i=0; i<n; i++) {
      printf("qxstart[%2d] = %23.16e\n", i, x[i]);
    }
#endif
//#define WRITE_QPS
#ifdef WRITE_QPS
    if (m0de==0) {
      FILE* fp = fopen("QPinit.dat", "w");
      fprintf(fp, "n = %d\n", n);
      fprintf(fp, "m = %d\n", m);
      fprintf(fp, "kmax = %d\n", kmax);
      fprintf(fp, "amax = %d\n" ,amax_);
      for (int i=1; i<=amax_; i++) {
        fprintf(fp, "a = %23.16e\n", a[i-1]);
      }
      int lamax = amax_ + m + 2;
      fprintf(fp, "lamax = %d\n" ,lamax);
      for (int i=1; i<=lamax+1; i++) {
        fprintf(fp, "la = %6d\n", la[i-1]);
      }
      for (int i=1; i<=n; i++) {
        fprintf(fp, "x = %23.16e\n", x[i-1]);
      }
      for (int i=1; i<=n+m; i++) {
        fprintf(fp, "bl = %23.16e\n", bl[i-1]);
      }
      for (int i=1; i<=n+m; i++) {
        fprintf(fp, "bu = %23.16e\n", bu[i-1]);
      }
      fprintf(fp, "fmin = %23.16e\n", fmin);
      fprintf(fp, "mlp = %6d\n", mlp);
      fprintf(fp, "mxws = %d\n", mxws);
      fprintf(fp, "mxlws = %d\n", mxlws);
      fclose(fp);
    }
    else {
      FILE* fp = fopen("QPbounds.dat", "w");
       fprintf(fp, "m0de = %d\n", m0de);
      for (int i=1; i<=n+m; i++) {
        fprintf(fp, "bl = %23.16e\n", bl[i-1]);
      }
      for (int i=1; i<=n+m; i++) {
        fprintf(fp, "bu = %23.16e\n", bu[i-1]);
      }
      fclose(fp);
    }
#endif
    F77_FUNC(bqpd,BQPD)(&n, &m, &k, &kmax, a, la, x, bl, bu, &f, &fmin,
        g, r, w, e, ls, alp, lp, &mlp, &peq, ws, lws,
        &m0de, &ifail, info, &iprint, &nout);
#ifdef TIME_BQPD
    times_.pivots += info[0];
#endif
    pivots_ += info[0];
    if(BqpdSolver::reinit_freq_ > 0 && haveHotStart_ && (ifail == 7 || ifail == 8) && m0de == 6){
      fprintf(stdout, "Reinitialize hot start...\n");
      copyFromHotStart();
      ifail = 0;
      F77_FUNC(bqpd,BQPD)(&n, &m, &k, &kmax, a, la, x, bl, bu, &f, &fmin,
			  g, r, w, e, ls, alp, lp, &mlp, &peq, ws, lws,
			  &m0de, &ifail, info, &iprint, &nout);
      printf("new ifail = %d\n", ifail);
    }
    if (ifail == 8 && BqpdSolver::m0de_ == 6) {
      fprintf(stdout, "Restarting Bqpd...");
      m0de = 0;
      tqp_->get_starting_point(n, 1, x, 0, NULL, NULL, m, 0, NULL);
      ifail = 0;
      F77_FUNC(bqpd,BQPD)(&n, &m, &k, &kmax, a, la, x, bl, bu, &f, &fmin,
			  g, r, w, e, ls, alp, lp, &mlp, &peq, ws, lws,
			  &m0de, &ifail, info, &iprint, &nout);
      printf("new ifail = %d\n", ifail);
    }
    while (ifail == 7) {
      // FIXME: For now, we just disable hot starts in case they were active
      if (haveHotStart_) unmarkHotStart();

      printf("Increasing fillin_factor from %e ", *fillin_factor_);
      *fillin_factor_ *= 2.;
      printf("to %e\n", *fillin_factor_);
      int mxws_new = kk + F77_FUNC(wsc,WSC).kkk + (5*n + (int)(*fillin_factor_*(double)amax_));
      real* ws_new = new real[mxws_new];
#ifdef InitializeAll
      for (int i=0; i<mxws_new; i++) {
	ws_new[i] = 42.;
      }
#endif
      CoinCopyN(ws, kk, ws_new);
      delete [] ws;
      ws = ws_new;
      mxws = mxws_new;
      F77_FUNC(wsc,WSC).mxws = mxws;

      m0de = 0;
      tqp_->get_starting_point(n, 1, x, 0, NULL, NULL, m, 0, NULL);
      ifail = 0;
      F77_FUNC(bqpd,BQPD)(&n, &m, &k, &kmax, a, la, x, bl, bu, &f, &fmin,
			  g, r, w, e, ls, alp, lp, &mlp, &peq, ws, lws,
			  &m0de, &ifail, info, &iprint, &nout);
    }
    if (ifail == 8) bad_warm_start_info_ = true;
#if 0
    for (int i=0; i<n; i++) {
      printf("qxsol[%2d] = %23.16e\n", i, x[i]);
    }
#endif
#if 0
    printf("ifail = %d\n", ifail);
    printf("final f = %e\n", f);
    printf("final f + obj_val = %e\n", f+tqp_->ObjVal());
#endif
#if 0
    int kkk = F77_FUNC(wsc,WSC).kkk;
    int lll = F77_FUNC(wsc,WSC).lll;
    printf("========= 3333333333 =============\n");
    printf("kk = %d ll = %d kkk = %d lll = %d mxws = %d mxlws = %d\n", kk, ll, kkk, lll, mxws, mxlws);
    for (int i=0; i<kk+kkk; i++) {
      printf("ws[%3d] = %15.8e\n ", i, ws[i]);
    }
    printf("--------\n");
    for (int i=kk+kkk; i<mxws; i++) {
      printf("ws[%3d] = %15.8e\n ", i, ws[i]);
    }
    for (int i=0; i<ll+lll; i++) {
      printf("lws[%5d] = %8d\n", i, lws[i]);
    }
    printf("------\n");
    for (int i=ll+lll; i<mxlws; i++) {
      printf("lws[%5d] = %8d\n", i, lws[i]);
    }
#endif
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
  {
    return errorNames_[0];
  }

  const std::string&
  BqpdSolver::UnsolvedBqpdError::solverName() const
  {
    return solverName_;
  }

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
    for (int i = 0 ; i < 14 ; i ++) {
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
  BqpdSolver::getEmptyWarmStart() const
  {
#if 0
    return new FilterWarmStart;
#endif
    printf("BqpdSolver::getEmptyWarmStart called \n");
    return NULL;
  }

  /** Check that warm start object is valid.*/
  bool 
  BqpdSolver::warmStartIsValid(const CoinWarmStart * ws) const{
    const BqpdWarmStart* bws = dynamic_cast<const BqpdWarmStart*>(ws);
    if (bws && ! bws->empty()) {
      return true;
    }
    return false;
  }

}//end namespace Bonmin

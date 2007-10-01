// (C) Copyright International Business Machines Corporation, Carnegie Mellon University 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 10/02/2006

#include "BonminConfig.h"

#include "BonFilterSolver.hpp"
#include "BonFilterWarmStart.hpp"

#include <fstream>

#include "CoinTime.hpp"
#include<algorithm>
typedef Bonmin::FilterSolver::fint fint;
typedef Bonmin::FilterSolver::real real;



typedef long ftnlen;
extern "C" {
void F77_FUNC(filtersqp,FILTERSQP)(
	fint *n, fint *m, fint *kmax, fint *maxa,
	fint *maxf, fint *mlp, fint *mxwk, fint *mxiwk,
	fint *iprint, fint *nout, fint *ifail, real *rho,
	real *x, real *c, real *f, real *fmin, real *bl,
	real *bu, real *s, real *a, fint *la, real *ws,
	fint *lws, real *lam, char *cstype, real *user,
	fint *iuser, fint *maxiter, fint *istat,
	real *rstat, ftnlen cstype_len);
}

//Static variables
static Ipopt::TNLP * tnlpSolved = NULL;
static fint nnz_h = -1;

static fint * hStruct = NULL;

//Permutation to apply to jacobian in order to get it row ordered
static int * permutationJac = NULL;
static int * permutationHess = NULL;
//static int * cache = NULL;


extern "C" {
//Access to filter common bloc
/* common block for problemname */
extern struct {
  fint  char_l;
  char pname[10];
} F77_FUNC(cpname,CPNAME);

/* common block for Hessian storage set to 0, i.e. NO Hessian */
extern struct {
  fint phl, phr, phc;
} F77_FUNC(hessc,HESSC);

/* common block for upper bound on filter */
extern struct {
  real ubd, tt;
} F77_FUNC(ubdc,UBDC);

/* common block for infinity & epslon */
extern struct {
  real infty, eps;
} F77_FUNC_(nlp_eps_inf,NLP_EPS_INF);

/* common block for prfinting from QP solver */
extern struct {
  fint n_bqpd_calls, n_bqpd_prfint;
} F77_FUNC_(bqpd_count,BQPD_COUNT);

/* common for scaling: scale_mode = 0 (none), 1 (variables), 2 (vars+cons) */
extern struct {
  fint scale_mode, phe;
} F77_FUNC(scalec,SCALEC);
}

extern "C" {

/// Objective function evaluation
void F77_FUNC(objfun,OBJFUN)(real *x, fint *n, real * f, real *user, fint * iuser, fint * errflag)
{
  (*errflag) = !tnlpSolved->eval_f(*n, x, 1, *f);
}

/** Constraint functions evaluation. */
void 
F77_FUNC(confun,CONFUN)(real * x, fint * n , fint *m, real *c, real *a, fint * la, real * user, fint * iuser,
	fint * errflag){
  (*errflag) = !tnlpSolved->eval_g(*n, x, 1, *m, c);
}

void
F77_FUNC(gradient,GRADIENT)(fint *n, fint *m, fint * mxa, real * x, real *a, fint * la,
	  fint * maxa, real * user, fint * iuser, fint * errflag){
  (*errflag) = !tnlpSolved->eval_grad_f(*n, x, 1, a);
  /// ATTENTION: Filter expect the jacobian to be ordered by row 
  int nnz = la[0] - *n - 1;
  double * values = new double [nnz];
  (*errflag) = !tnlpSolved->eval_jac_g(*n, x, 1, *m, nnz, NULL, NULL, values) || (*errflag);
  a+= *n;
  for(int i = 0 ; i < nnz ; i++)
    {
      int indice = permutationJac[i];
      if(indice > nnz)
	{
#ifndef NDEBUG
	  std::cout<<"Error in gradient computation, i: "<<i
                   <<" in row order "<<permutationJac[i]<<std::endl;
#endif
         }
      *a++ = values[indice];
    }
  delete [] values;
}

/* evaluation of the Hessian of the Lagrangian */
 void
F77_FUNC(hessian,HESSIAN)(real *x, fint *n, fint *m, fint *phase, real *lam,
	     real *ws, fint *lws, real *user, fint *iuser,
	     fint *l_hess, fint *li_hess, fint *errflag)
{
  Number obj_factor = (*phase == 1)? 0. : 1.;
  fint  end = nnz_h + (*n)  + 2;

  for(int i = 0 ; i < end ; i++)
    {
      lws[i] = hStruct[i];
    }
  *l_hess = nnz_h;
  *li_hess = nnz_h + *n + 3;
  Number * mlam = NULL;
  if (*m > 0) {
    mlam = new Number[*m];
  }
  for (int i = 0; i<*m; i++) {
    mlam[i] = -lam[*n+i];
  }
  Number * values = new Number [nnz_h];
  (*errflag) = !tnlpSolved->eval_h(*n, x, 1, obj_factor, *m, mlam ,1, hStruct[0] - 1, NULL, NULL, values);
  delete [] mlam;
  for(int i = 0 ; i < nnz_h ; i++) ws[i] = values[permutationHess[i]];
  delete [] values;
}

}

namespace Bonmin{

  struct Transposer{
    const Index* rowIndices;
    const Index* colIndices;
    bool operator()(int i, int j){
      return rowIndices[i]<rowIndices[j] ||
	(rowIndices[i]==rowIndices[j] && colIndices[i] < colIndices[j]);}
  };

  // Convert a sparse matrix from triplet format to row ordered packed matrix
  void FilterSolver::TMat2RowPMat(bool symmetric, fint n, fint m, int nnz,
                                  const Index* iRow,
				  const Index* iCol, int * permutation2,
				  fint * lws, int nnz_offset, int n_offset,
                                  Ipopt::TNLP::IndexStyleEnum index_style)
  {
    for(int i = 0 ; i < nnz ; i++)
      permutation2[i] = i;

    Transposer lt;
    if (symmetric) {
      Index* tmpRow = new Index[nnz];
      Index* tmpCol = new Index[nnz];
      for (int i=0; i<nnz; i++) {
        const Index& irow = iRow[i];
        const Index& jcol = iCol[i];
        if (irow > jcol) {
          tmpRow[i] = irow;
          tmpCol[i] = jcol;
        }
        else {
          tmpRow[i] = jcol;
          tmpCol[i] = irow;
        }
      }
      lt.rowIndices = tmpRow;
      lt.colIndices = tmpCol;
    }
    else {
      lt.rowIndices = iRow;
      lt.colIndices = iCol;
    }

    std::sort(permutation2, &permutation2[nnz], lt);

    const int idx_offset = (index_style == Ipopt::TNLP::C_STYLE);
    fint row = 1-idx_offset;
    lws[0] = nnz + nnz_offset + 1;
    fint * inds = lws + nnz_offset + 1;
    fint * start = inds + nnz + n_offset;
    *start++ = 1 + nnz_offset;
    for(fint i = 0 ; i < nnz ; i++)
      {
	inds[i] = lt.colIndices[permutation2[i]] + idx_offset;
	//DBG_ASSERT(RowJac[permutation2[i]] >= row);
	if(lt.rowIndices[permutation2[i]] > row) {
	  for(;row < lt.rowIndices[permutation2[i]] ; row++)
	    *start++ = i + nnz_offset + 1;
	}
      }
    for(;row <= m-idx_offset ; row++)
      *start++ = nnz + nnz_offset +1;

#if 0
    for (int i = 0; i<nnz_offset+1; i++)
      printf("lws[%3d] = %3d\n", i, lws[i]);
    for (int i = nnz_offset+1; i<nnz_offset+nnz+1; i++)
      printf("lws[%3d] = %3d  [%3d,%3d]\n", i, lws[i], lt.rowIndices[permutation2[i-nnz_offset-1]], lt.colIndices[permutation2[i-nnz_offset-1]]);
    for (int i = nnz_offset+nnz+1; i<lws[0]+m+2; i++)
      printf("lws[%3d] = %3d\n", i, lws[i]);
#endif

    if (symmetric) {
      delete [] lt.rowIndices;
      delete [] lt.colIndices;
    }


  }

#if 0
  // Convert a sparse matrix from triplet format to row ordered packed matrix
  void TMat2RowPMat_(fint n, fint m, int nnz, const Index* iRow,
				  const Index* iCol, int * permutation2,
				  fint * lws, int offset)
  {
    for(int i = 0 ; i < nnz ; i++)
      permutation2[i] = i;

    Transposer lt;
    lt.rowIndices = iRow;
    lt.colIndices = iCol;

    std::sort(permutation2, &permutation2[nnz], lt);

    fint row = 1;
    lws[0] = nnz + offset + 1;
    fint * inds = lws + offset + 1;
    fint * start = inds + nnz + 1;

    for(fint i = 0 ; i < nnz ; i++)
      {
	inds[i] = iCol[permutation2[i]];
	//DBG_ASSERT(RowJac[permutation2[i]] >= row);
	if(iRow[permutation2[i]] >= row) {
	  for(;row <= iRow[permutation2[i]] ; row++)
	    *start++ = i + offset + 1;
	}
      }
    for(;row <= m+1 ; row++)
      *start++ = nnz + offset +1;

    //deleteme
    for (int i = 0; i<offset+1; i++)
      printf("alws[%3d] = %3d\n", i, lws[i]);
    for (int i = offset+1; i<offset+nnz+1; i++)
      printf("alws[%3d] = %3d  [%3d,%3d]\n", i, lws[i], lt.rowIndices[permutation2[i-offset-1]], lt.colIndices[permutation2[i-offset-1]]);
    for (int i = offset+nnz+1; i<lws[0]+m+2; i++)
      printf("alws[%3d] = %3d\n", i, lws[i]);
  }
  

  // Convert a sparse matrix from triplet format to row ordered packed matrix
  void TMat2ColPMat_(fint n, fint m, int nnz, const Index* iRow,
				  const Index* iCol,
				  int* permutationHess2, fint * lws, int offset)
  {
    for(int i = 0 ; i < nnz ; i++)
      permutationHess2[i] = i;

    fint col = 1;
    lws[0] = nnz + offset + 1;
    fint * inds = lws + 1;
    fint * start = inds + nnz + offset;
    
    Transposer lt;
    lt.rowIndices = iCol;
    lt.colIndices = iRow;

    std::sort(permutationHess2, permutationHess2 + nnz, lt);

   for(fint i = 0 ; i < nnz ; i++)
      {
	inds[offset + i] = iRow[permutationHess2[i]];
	if(iCol[permutationHess2[i]] >= col) {
	  for(;col <= iCol[permutationHess2[i]] ; col++)
	    *start++ = i + offset + 1;
	}
      }
    for(;col <= n+1 ; col++)
      *start++ = nnz + offset +1;

    //deleteme
    for (int i = 0; i<offset+1; i++)
      printf("blws[%3d] = %3d\n", i, lws[i]);
    for (int i = offset+1; i<offset+nnz+1; i++)
      printf("blws[%3d] = %3d  [%3d,%3d]\n", i, lws[i], lt.rowIndices[permutationHess2[i-offset-1]], lt.colIndices[permutationHess2[i-offset-1]]);
    for (int i = offset+nnz+1; i<lws[0]+m+2; i++)
      printf("blws[%3d] = %3d\n", i, lws[i]);
  }
#endif


  std::string FilterSolver::solverName_ = "filter SQP";

void
FilterSolver::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
  roptions->SetRegisteringCategory("FilterSQP options", RegisteredOptions::FilterCategory);
  roptions->AddLowerBoundedNumberOption("eps", "Tolerance for SQP solver",
					0., 1, 1e-06, "");
  roptions->AddLowerBoundedNumberOption("infty","A large number (1E20)",0.,1, 1e20, "");
  roptions->AddBoundedIntegerOption("iprint", "Print level (0=silent, 3=verbose)", 0,6,0);
  roptions->AddLowerBoundedIntegerOption("kmax", "Dimension of null-space",
				 -1, -1, "");
  roptions->AddLowerBoundedIntegerOption("maxf","Maximum filter length",0,100);
  roptions->AddLowerBoundedIntegerOption("maxiter", "Maximum number of iterations",0,1000);
  roptions->AddLowerBoundedIntegerOption("mlp","Maximum level for degeneracy (bqpd)",0, 1000);
  roptions->AddLowerBoundedIntegerOption("mxlws", "FINTEGER workspace increment", 0, 500000);
  roptions->AddLowerBoundedIntegerOption("mxws", "REAL workspace increment",
					0,2000000);
  roptions->AddLowerBoundedNumberOption("rho", "Initial trust region size",0,1,10.);
  //  roption->AddLowerBoundedIntegerOption("timing", "whether to time evaluations (1 = yes)");
  roptions->AddLowerBoundedNumberOption("tt", "Parameter for upper bound on filter",0,1, 125e-2);
  roptions->AddLowerBoundedNumberOption("ubd", "Parameter for upper bound on filter", 0 , 1,1e2);

 }


  FilterSolver::FilterSolver(bool createEmpty /* = false */)
    :
    TNLPSolver(),
    cached_(NULL)
{
}

FilterSolver::FilterSolver(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
                           Ipopt::SmartPtr<Ipopt::OptionsList> options,
                           Ipopt::SmartPtr<Ipopt::Journalist> journalist):
    TNLPSolver(roptions, options, journalist),
    cached_(NULL)
{
}

                           
Ipopt::SmartPtr <TNLPSolver>
FilterSolver::clone(){
  Ipopt::SmartPtr<FilterSolver> retval = new FilterSolver(true);
  *retval->options_ = *options_; // Copy the options
  retval->roptions_ = roptions_; // only copy pointers of registered options
  retval->journalist_ = journalist_; // and journalist
  return GetRawPtr(retval);
}

FilterSolver::~FilterSolver(){
}

bool
FilterSolver::Initialize(std::string optFile){
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
FilterSolver::Initialize(std::istream &is){

  Index ivalue;
  options_->GetIntegerValue("print_level", ivalue, "");
  EJournalLevel print_level = (EJournalLevel)ivalue;
  SmartPtr<Journal> stdout_jrnl = journalist_->GetJournal("console");
  if (IsValid(stdout_jrnl)) {
    // Set printlevel for stdout
    stdout_jrnl->SetAllPrintLevels(print_level);
    stdout_jrnl->SetPrintLevel(J_DBG, J_NONE);
  }

  if(is.good()){
    options_->ReadFromStream(*journalist_, is);
  }
  return true;
}

/// Solves a problem expresses as a TNLP 
TNLPSolver::ReturnStatus 
FilterSolver::OptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> & tnlp)
{
  if (IsNull(cached_) || !cached_->use_warm_start_in_cache_) {
    cached_ = new cachedInfo(tnlp, options_);
  }
  return callOptimizer();
}

/// Solves a problem expresses as a TNLP 
TNLPSolver::ReturnStatus 
FilterSolver::ReOptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> & tnlp)
{
  assert(tnlp == cached_->tnlp_);
  //rescan bounds which may have changed
  assert(cached_->bounds);
  int n = cached_->n;
  int m = cached_->m;
  tnlp->get_bounds_info(n, cached_->bounds, &cached_->bounds[n+m], 
			m, &cached_->bounds[n], &cached_->bounds[2*n + m]);

  tnlpSolved = static_cast<Ipopt::TNLP *>(Ipopt::GetRawPtr(tnlp));
  nnz_h = cached_->nnz_h_;

  hStruct = cached_->hStruct_;

//Permutation to apply to jacobian in order to get it row ordered
  permutationJac = cached_->permutationJac_;
  permutationHess = cached_->permutationHess_;


  return callOptimizer();
}



void
FilterSolver::cachedInfo::initialize(const Ipopt::SmartPtr<Ipopt::TNLP> & tnlp,
				     Ipopt::SmartPtr<Ipopt::OptionsList>& options)
{
  // 1) Get some dimensions
  // 1.a) First from ampl
  int  nnz_jac_g;

  Ipopt::TNLP::IndexStyleEnum index_style;
  Ipopt::Index nv, nc, nnz_j, nnz_hess;
  tnlp->get_nlp_info( nv, nc, 
                     nnz_j, (Ipopt::Index&) nnz_hess, 
                     index_style);
  n = nv;
  m = nc;
  nnz_jac_g = nnz_j;
  nnz_h_ = nnz_hess;

  nnz_h = nnz_h_;
  

  // 1.b) then from options
  Ipopt::Index kmax_ipt;
  options->GetIntegerValue("kmax", kmax_ipt, "filter.");
  if (kmax_ipt == -1) {
    kmax = n;
  }
  else {
    kmax = kmax_ipt;
    kmax = min(kmax,n);
  }
  Ipopt::Index mlp_ipt;
  options->GetIntegerValue("mlp", mlp_ipt,"filter.");
  mlp = mlp_ipt;

  Ipopt::Index  maxf_ipt;
  options->GetIntegerValue("maxf", maxf_ipt,"filter.");
  maxf = maxf_ipt;

  fint mxwk0;
  Ipopt::Index mxwk0_ipt;
  options->GetIntegerValue("mxws", mxwk0_ipt, "filter.");
  mxwk0 = mxwk0_ipt;

  fint mxiwk0;
  Ipopt::Index mxiwk0_ipt;
  options->GetIntegerValue("mxlws",  mxiwk0_ipt, "filter.");
  mxiwk0 = mxiwk0_ipt;
  // Setup storage for Filter
  int nplusm = n + m;
  //Starting point
  x = new real [n];

  //tnlp->get_starting_point(n, 1, x, 0, NULL, NULL, m, 0, NULL);
  use_warm_start_in_cache_ = false;
  //for(int i = 0 ; i < n ; i++) x[i] = 0;
  lam = new real [n+m];
  //#define InitializeAll
#ifdef InitializeAll
  for(int i = 0 ; i < n+m ; i++) lam[i] = 0.; 
#endif
  //bounds
  bounds = new real [2*n + 2*m];
  
  tnlp->get_bounds_info(n, bounds, bounds + nplusm, m, bounds + n, bounds + n + nplusm);

#if 0
  double infty = F77_FUNC_(nlp_eps_inf,NLP_EPS_INF).infty;
  // AW: I don't think we need this, it isn't done for ReOptimize either
  for(int i = 0 ; i < nplusm ; i++){
  if(bounds[i] < -infty) bounds[i] = - infty;}

  real * ubounds = bounds + nplusm;
  for(int i = 0 ; i < nplusm ; i++){
  if(ubounds[i] > infty) ubounds[i] = infty;}
#endif

  maxa = n + nnz_jac_g;
  fint maxia = n + nnz_jac_g + m + 3;
  a = new real[maxa];
  la = new fint [maxia];

  int * RowJac = new int [nnz_jac_g];  
  int * ColJac = new int [nnz_jac_g];
  
  la[nnz_jac_g + n + 1] = 1;

  for(fint i = 1; i <= n ; i++)
    la[i] = i;// - (index_style == Ipopt::TNLP::C_STYLE);
  tnlp->eval_jac_g(  nv, NULL, 0, nc , nnz_j,  RowJac,  ColJac, NULL);

  permutationJac = permutationJac_ = new int [nnz_jac_g];
  TMat2RowPMat(false, n, m, nnz_jac_g,  RowJac, ColJac, permutationJac,
	       la, n, 1, index_style);

  delete [] RowJac;
  delete [] ColJac;

  // Now setup hessian
  permutationHess = permutationHess_ = new int[nnz_h];
  hStruct_ = new fint[nnz_h + n + 3];
  int * cache = new int[2*nnz_h + 1];
  F77_FUNC(hessc,HESSC).phl = 1;
  tnlp->eval_h((Ipopt::Index&) n, NULL, 0, 1., (Ipopt::Index&) m, NULL, 0, (Ipopt::Index&) nnz_h, cache + nnz_h, cache  , NULL);

  TMat2RowPMat(true, n, n, nnz_h, cache, cache + nnz_h, permutationHess,
	       hStruct_, 0, 0, index_style);

  delete [] cache;
  // work arrays
  fint lh1 = nnz_h + 8 + 2 * n + m;
  maxWk = 21*n + 8*m + mlp + 8*maxf + lh1 + kmax*(kmax+9)/2 + mxwk0;
  maxiWk = 13*n + 4*m + mlp + lh1 + kmax + 113 + mxiwk0;

  ws = new real[maxWk];
  lws = new fint[maxiWk];
#ifdef InitializeAll
  for(int i = 0 ; i < maxWk ; i++) ws[i] = 0;
  for(int i = 0 ; i < maxiWk ; i++) lws[i] = 0;
#endif

  // Setup global variables and static variables
  hStruct = hStruct_;
  tnlpSolved = static_cast<Ipopt::TNLP *>(Ipopt::GetRawPtr(tnlp));

  options->GetNumericValue("ubd",F77_FUNC(ubdc,UBDC).ubd, "filter.");
  options->GetNumericValue("tt", F77_FUNC(ubdc,UBDC).tt, "filter.");
  options->GetNumericValue("eps", F77_FUNC_(nlp_eps_inf,NLP_EPS_INF).eps, "filter.");
  options->GetNumericValue("infty", F77_FUNC_(nlp_eps_inf,NLP_EPS_INF).infty, "filter.");
  rho = 10.;
  maxiter = 1000;
  options->GetIntegerValue("maxiter", (Ipopt::Index &) maxiter, "filter.");
  options->GetNumericValue("rho",rho,"filter.");


  // Set up scaling
  F77_FUNC(scalec,SCALEC).scale_mode = 0;
  s = new real [n+m];
  
  istat = new fint[14];
  rstat = new real[7];
#ifdef InitializeAll
  for (int i=0; i<14; i++) {
    istat[0] = 43;
  }
  for (int i=0; i<7; i++) {
    rstat[0] = 42.;
  }
#endif

  fmin = -1e100;
  Ipopt::Index bufy;
  options->GetIntegerValue("iprint",bufy, "filter.");
  iprint = bufy;
  nout = 6;
  cstype = new char[m];
  Ipopt::TNLP::LinearityType * const_types =
        new Ipopt::TNLP::LinearityType[m];
  tnlp->get_constraints_linearity(m, const_types);
  for(int i = 0 ; i < m ; i++)
  {
    if(const_types[i] == Ipopt::TNLP::LINEAR){
     cstype[i] = 'L';}
    else 
     cstype[i] = 'N';
  } 
  delete [] const_types;
  c = new double[m];
  tnlp_ = Ipopt::GetRawPtr(tnlp);
}

/// Solves a problem expresses as a TNLP 
TNLPSolver::ReturnStatus 
FilterSolver::callOptimizer()
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
  case 4:
    optimizationStatus = TNLPSolver::provenInfeasible;
    status = Ipopt::LOCAL_INFEASIBILITY;
    break;
  case 5:
    optimizationStatus = TNLPSolver::computationError;
    status = Ipopt::INTERNAL_ERROR;
    break;
  case 6:
    optimizationStatus = TNLPSolver::iterationLimit;
    status = Ipopt::MAXITER_EXCEEDED;
    break;
  case 7:
    optimizationStatus = TNLPSolver::externalException;
    status = Ipopt::INTERNAL_ERROR;
    break;
  case 8:
  case 9:
  case 10:
    optimizationStatus = TNLPSolver::exception;
    status = Ipopt::INTERNAL_ERROR;
    break;
  }

  Number* mlam = NULL;
  if (cached_->m>0) {
    mlam = new Number[cached_->m];
  }
  for (int i = 0; i<cached_->m; i++) {
    mlam[i] = -cached_->lam[cached_->n+i];
  }
  Number* z_L = new Number[cached_->n];
  Number* z_U = new Number[cached_->n];
  const int os = cached_->n+cached_->m;
  for (int i=0; i<cached_->n; i++) {
    if (cached_->x[i] == cached_->bounds[i]) {
      z_L[i] = Max(0.,cached_->lam[i]);
    }
    else {
      z_L[i] = 0.;
    }
    if (cached_->x[i] == cached_->bounds[os+i]) {
      z_U[i] = Max(0.,-cached_->lam[i]);
    }
    else {
      z_U[i] = 0.;
    }
  }
  cached_->tnlp_->finalize_solution(status, cached_->n, 
			  cached_->x, z_L, z_U, 
			  cached_->m, cached_->c, mlam, 
			  cached_->f, NULL, NULL);
  delete [] mlam;
  delete [] z_L;
  delete [] z_U;
  return optimizationStatus;
}

/** Optimize problem described by cache with filter.*/
void 
FilterSolver::cachedInfo::optimize()
{
  if (use_warm_start_in_cache_) {
    ifail = -1;
    use_warm_start_in_cache_ = false;
  }
  else {
    tnlp_->get_starting_point(n, 1, x, 0, NULL, NULL, m, 0, NULL);
    ifail = 0;
  }
  printf("ifail = %d\n", ifail);
  cpuTime_ = - CoinCpuTime();
  fint cstype_len = 1;
  rho = 10;
  //  rho = 1e6;
  //  printf("rho = %e\n", rho);
#if 0
  printf("========= 3333333333333333 =============\n");
  for (int i=0; i<n; i++) {
    printf("xL[%3d] = %15.8e  xU[%3d] = %15.8e\n", i, bounds[i], i, bounds[m+n+i]);
  }
#endif
  F77_FUNC(filtersqp,FILTERSQP)(&n, &m, &kmax, & maxa, &maxf, &mlp, &maxWk, 
	     &maxiWk, &iprint, &nout, &ifail, &rho, x, 
	     c, &f, &fmin, bounds, 
	     bounds + n + m, 
	     s, a, la, 
	     ws, lws, lam, cstype, 
	     NULL, NULL,
	     &maxiter, istat, rstat,
	     cstype_len);
   cpuTime_ += CoinCpuTime();
}

std::string
FilterSolver::UnsolvedFilterError::errorNames_[1] =
{"Internal error in Filter SQP."};

std::string
FilterSolver::UnsolvedFilterError::solverName_ = 
          "filterSqp";

const std::string& 
FilterSolver::UnsolvedFilterError::errorName() const
{ return errorNames_[0];}

const std::string& 
FilterSolver::UnsolvedFilterError::solverName() const
{ return solverName_;}

bool
FilterSolver::setWarmStart(const CoinWarmStart * warm, 
			   Ipopt::SmartPtr<TMINLP2TNLP> tnlp)
{
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
  return true;
}

CoinWarmStart *
FilterSolver::getWarmStart(Ipopt::SmartPtr<TMINLP2TNLP> tnlp) const
{
  return new FilterWarmStart(cached_->n, cached_->x,
			     cached_->n+cached_->m, cached_->lam,
			     cached_->maxiWk, cached_->lws, cached_->istat);
}

CoinWarmStart * 
FilterSolver::getEmptyWarmStart() const {
  return new FilterWarmStart;}



}//end namespace Bonmin

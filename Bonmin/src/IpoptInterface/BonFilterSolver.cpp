#include "BonFilterSolver.hpp"
#include <fstream>



typedef Bonmin::FilterSolver::fint fint;
typedef Bonmin::FilterSolver::real real;



typedef long ftnlen;

void filtersqp_(
	fint *n, fint *m, fint *kmax, fint *maxa,
	fint *maxf, fint *mlp, fint *mxwk, fint *mxiwk,
	fint *iprint, fint *nout, fint *ifail, real *rho,
	real *x, real *c, real *f, real *fmin, real *bl,
	real *bu, real *s, real *a, fint *la, real *ws,
	fint *lws, real *lam, char *cstype, real *user,
	fint *iuser, fint *maxiter, fint *istat,
	real *rstat, ftnlen cstype_len);


//Static variables
static Ipopt::TNLP * tnlpSolved = NULL;
static fint nnz_h = -1;

static fint * hStruct = NULL;
//static int * cache = NULL;


//Access to filter common bloc
/* common block for problemname */
extern struct {
  fint  char_l;
  char pname[10];
} cpname_;

/* common block for Hessian storage set to 0, i.e. NO Hessian */
extern struct {
  fint phl, phr, phc;
} hessc_;

/* common block for upper bound on filter */
extern struct {
  real ubd, tt;
} ubdc_;

/* common block for infinity & epslon */
extern struct {
  real infty, eps;
} nlp_eps_inf__;

/* common block for prfinting from QP solver */
extern struct {
  fint n_bqpd_calls, n_bqpd_prfint;
} bqpd_count__;

/* common for scaling: scale_mode = 0 (none), 1 (variables), 2 (vars+cons) */
extern struct {
  fint scale_mode, phe;
} scalec_;


extern "C" {

/// Objective function evaluation
void objfun_(real *x, fint *n, real * f, real *user, fint * iuser, fint * errflag)
{
  tnlpSolved->eval_f(*n, x, 1, *f);
}

/** Constraint functions evaluation. */
void 
confun_(real * x, fint * n , fint *m, real *c, real *a, fint * la, real * user, fint * iuser,
	fint * errflags){
  tnlpSolved->eval_g(*n, x, 1, *m, c);
}

void
gradient_(fint *n, fint *m, fint * mxa, real * x, real *a, fint * la,
	  fint * maxa, real * user, fint * iuser, fint * errflag){
  tnlpSolved->eval_grad_f(*n, x, 1, a);
  /// ATTENTION: Filter expect the jacobian to be ordered by row 
  tnlpSolved->eval_jac_g(*n, x, 1, *m, la[0] - *n, NULL, NULL, &a[*n]);
}

/* evaluation of the Hessian of the Lagrangian */
 void
hessian_(real *x, fint *n, fint *m, fint *phase, real *lam,
	     real *ws, fint *lws, real *user, fint *iuser,
	     fint *l_hess, fint *li_hess, fint *errflag)
{
  real obj_factor = (*phase == 1)? 1. : 0.;
  fint * begin = hStruct;
  fint * end = begin + 2 * begin[0];
  *lws++ = *begin++;

  for( ; begin != end ; lws += 2, begin += 2)
    {
      lws[0] = begin[0];
      lws[1] = begin[1];
    }
  tnlpSolved->eval_h(*n, x, 1, obj_factor, *m, &lam[*n],1, hStruct[0], NULL, NULL, ws);
}

}

namespace Bonmin{

void
FilterSolver::registerAllOptions(){
  roptions_->SetRegisteringCategory("FilterSQP options");
  roptions_->AddLowerBoundedNumberOption("eps", "Tolerance for SQP solver",
					0., 1, 1e-06, "");
  roptions_->AddLowerBoundedNumberOption("infty","A large number (1E20)",0.,1, 1e20, "");
  roptions_->AddBoundedIntegerOption("iprint", "Print level (0=silent, 3=verbose)", 0,3,0);
  roptions_->AddLowerBoundedIntegerOption("kmax", "Dimension of null-space",
					 0,500, "");
  roptions_->AddLowerBoundedIntegerOption("maxf","Maximum filter length",0,100);
  roptions_->AddLowerBoundedIntegerOption("maxiter", "Maximum number of iterations",0,1000);
  roptions_->AddLowerBoundedIntegerOption("mlp","Maximum level for degeneracy (bqpd)",0, 1000);
  roptions_->AddLowerBoundedIntegerOption("mxlws", "FINTEGER workspace increment", 0, 500000);
  roptions_->AddLowerBoundedIntegerOption("mxws", "REAL workspace increment",
					0,200000);
  roptions_->AddLowerBoundedNumberOption("rho", "Initial trust region size",0,1,10.);
  //  roption->AddLowerBoundedIntegerOption("timing", "whether to time evaluations (1 = yes)");
  roptions_->AddLowerBoundedNumberOption("tt", "Parameter for upper bound on filter",0,1, 125e-2);
  roptions_->AddLowerBoundedNumberOption("ubd", "Parameter for upper bound on filter", 0 , 1,1e2);
}


FilterSolver::FilterSolver():
  journalist_(new Ipopt::Journalist()),
  options_(new Ipopt::OptionsList()),
  roptions_(new Ipopt::RegisteredOptions()),
  cached_(NULL)
{
  try{
  Ipopt::SmartPtr<Ipopt::Journal> stdout_journal =
    journalist_->AddFileJournal("console", "stdout", Ipopt::J_ITERSUMMARY);
  
  registerAllOptions();

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
    Ipopt::IpoptException E("Uncaught exception in FilterSolver::FilterSolver()",
		      "BonFilterSolver.cpp",-1);
    throw E;
  }
}

FilterSolver::~FilterSolver(){
  cached_ = NULL;
}

void 
FilterSolver::Initialize(std::string optFile){
    std::ifstream is;
    if (optFile != "") {
      try {
        is.open(optFile.c_str());
      }
      catch(std::bad_alloc) {
        journalist_->Printf(Ipopt::J_SUMMARY, Ipopt::J_MAIN, "\nEXIT: Not enough memory.\n");
        throw -1;
      }
      catch(...) {
        Ipopt::IpoptException E("Unknown Exception caught in ipopt", "Unknown File", -1);
        E.ReportException(*journalist_);
        throw -1;
      }
    }
    Initialize(is);
    if (is) {
      is.close();
    }
  
}

void
FilterSolver::Initialize(std::istream &is){
  if(is.good()){
    options_->ReadFromStream(*journalist_, is);
  }
}

/// Solves a problem expresses as a TNLP 
TNLPSolver::ReturnStatus 
FilterSolver::OptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> & tnlp)
{
  cached_ = new cachedInfo(tnlp, options_);
  return callOptimizer();
}

/// Solves a problem expresses as a TNLP 
TNLPSolver::ReturnStatus 
FilterSolver::ReOptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> & tnlp)
{
  assert(tnlp == cached_->tnlp_);
  cached_->ifail = -1;
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
  tnlp->get_nlp_info((Ipopt::Index&) n,(Ipopt::Index&) m, 
                     (Ipopt::Index&) nnz_jac_g, (Ipopt::Index&) nnz_h, 
                     index_style);


  

  // 1.b) then from options
  options->GetIntegerValue("kmax", (Ipopt::Index&) kmax, "filter.");

  options->GetIntegerValue("mlp", (Ipopt::Index&) mlp,"filter.");

  options->GetIntegerValue("maxf",(Ipopt::Index&) maxf,"filter.");

  fint mxwk0;
  options->GetIntegerValue("mxws", (Ipopt::Index&) mxwk0, "filter.");

  fint mxiwk0;
  options->GetIntegerValue("mxlws", (Ipopt::Index&) mxiwk0, "filter.");

  // Setup storage for Filter

  //Starting point
  x = new real [n];

  tnlp->get_starting_point(n, 1, x, 0, NULL, NULL, m, 0, NULL);
  
  lam = new real [n+m];
  
  //bounds
  bounds = new real [2*n + 2*m];
  
  tnlp->get_bounds_info(n, bounds, &bounds[n+m], m, &bounds[n], &bounds[2*n + m]);
  maxa = n + nnz_jac_g;
  fint maxia = n + nnz_jac_g + m + 3;
  a = new real[maxa];
  la = new fint [maxia];
  fint * start = & la[nnz_jac_g + 1];
  int * RowJac = new int [nnz_jac_g];  
  int * ColJac = new int [nnz_jac_g];
  
  la[0] = nnz_jac_g + n + 1;

  *start++ = 1;

  for(fint i = 1; i <= n ; i++)
    la[i] = i;// - (index_style == Ipopt::TNLP::C_STYLE);
  tnlp->eval_jac_g( (int) n, NULL, 0,(int) m , (int) nnz_jac_g,  RowJac,  ColJac, NULL);

  fint row = -1;
  for(fint i = 0 ; i < nnz_jac_g ; i++)
    {
      la[1+n+i] = ColJac[i];
      DBG_ASSERT(RowJac[i] >= row);
      if(RowJac[i] > row) {
	for(row++; row <= RowJac[i] ; row++)
	  *start++ = i + n + 1;
      }
    }
  for(row++; row <= m ; row++)
    *start++ = nnz_jac_g + n +1;

  delete [] RowJac;
  delete [] ColJac;

  // Now setup hession
  
  hStruct_ = new fint[2 * nnz_h + 1];
  int * cache = new int[2*nnz_h];
  hStruct_[0] = nnz_h + 1;
  hessc_.phl = 1;
  tnlp->eval_h((Ipopt::Index&) n, NULL, 0, 1., (Ipopt::Index&) m, NULL, 0, (Ipopt::Index&) nnz_h, cache , cache + nnz_h , NULL);

  int end = 2 * nnz_h;
  for(int i = 1 ; i <= end ; i++)
    {
      hStruct_[i] = cache[i];
    }

  delete [] cache;
  // work arrays
  fint lh1 = nnz_h + 8 + 2 * n + m;
  maxWk = 21*n + 8*m + mlp + 8*maxf + lh1 + kmax*(kmax+9)/2 + mxwk0;
  maxiWk = 13*n + 4*m + mlp + lh1 + kmax + 113 + mxiwk0;

  ws = new real[maxWk];
  lws = new fint[maxiWk];

  // Setup global variables and static variables
  hStruct = hStruct_;
  tnlpSolved = static_cast<Ipopt::TNLP *>(Ipopt::GetRawPtr(tnlp));

  options->GetNumericValue("ubd",ubdc_.ubd, "filter.");
  options->GetNumericValue("tt", ubdc_.tt, "filter.");
  options->GetNumericValue("eps", nlp_eps_inf__.eps, "filter.");
  options->GetNumericValue("infty", nlp_eps_inf__.infty, "filter.");
  rho = 10.;
  maxiter = 1000;
  options->GetIntegerValue("maxiter", (Ipopt::Index &) maxiter, "filter.");
  options->GetNumericValue("rho",rho,"filter.");


  istat = new fint[14];
  rstat = new real[7];

  fmin = -1e100;
  options->GetIntegerValue("iprint",(Ipopt::Index &) iprint, "filter.");
  cstype = new char[m];
  for(int i = 0 ; i < m ; i++)
     cstype[i] = 'N'; 
  c = new double[m];
  tnlp_ = Ipopt::GetRawPtr(tnlp);
}

/// Solves a problem expresses as a TNLP 
TNLPSolver::ReturnStatus 
FilterSolver::callOptimizer()
{
  cached_->optimize();
  TNLPSolver::ReturnStatus optimizationStatus;
  Ipopt::SolverReturn status;
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

  cached_->tnlp_->finalize_solution(status, cached_->n, 
			  cached_->x, cached_->lam, cached_->lam, 
			  cached_->m, cached_->c, cached_->lam + cached_->n, 
			  cached_->f);
  return optimizationStatus;
}

/** Optimize problem described by cache with filter.*/
void 
FilterSolver::cachedInfo::optimize()
{
  fint cstype_len = 1;
  filtersqp_(&n, &m, &kmax, & maxa, &maxf, &mlp, &maxWk, 
	     &maxiWk, &iprint, &nout, &ifail, &rho, x, 
	     c, &f, &fmin, bounds, 
	     bounds + n + m, 
	     NULL, a, la, 
	     ws, lws, lam, cstype, 
	     NULL, NULL,
	     &maxiter, istat, rstat,
	     cstype_len);
}


}//end namespace Bonmin

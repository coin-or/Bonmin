// (C) Copyright International Business Machines Corporation, 2006
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 10/02/2006

//Bonmin includes
#include "BonTNLPSolver.hpp"

//Ipopt includes
#include "IpBlas.hpp"

//Standard includes
#include <fstream>

#include "BonColReader.hpp"


namespace Bonmin{
  using namespace Ipopt;

  
#if 0
  string TNLPSolver::returnCodes[TNLPSolver::numReturnCodes] = {
    "Hit iteration limit" /**iterationLimit*/,
    "Some error was made in computations." /**computationError*/,
    "Problem has more equations than free variables." /** notEnoughFreedom*/,
    "Problem is not well defined."/**illDefinedProblem*/,
    "Illegal option set."/**illegalOption*/,
    "Exception in some the third-party code used by solver." /**externalException*/,
    "Unrecognized exception."/**exception*/,
    "Problem solved to optimality"/**solvedOptimal*/,
    "Problem solvedd to acceptable level of tolerance"/**solvedOptimalTol*/,
    "Problem is infeasible" /**provenInfeasible*/,
    "Problem is unbounded", /**unbounded*/};
#endif
    
    
  TNLPSolver::TNLPSolver():
    start_time_(0),
    time_limit_(DBL_MAX)
  {
    initializeOptionsAndJournalist();
  }

  TNLPSolver::TNLPSolver(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
    Ipopt::SmartPtr<Ipopt::OptionsList> options,
    Ipopt::SmartPtr<Ipopt::Journalist> journalist,
    const std::string & prefix):
    journalist_(journalist),
    options_(options),
    roptions_(roptions),
    prefix_(prefix),
    start_time_(0),
    time_limit_(DBL_MAX)
  {
  }

  TNLPSolver::TNLPSolver(const TNLPSolver &other):
    journalist_(other.journalist_),
    options_(NULL),
    roptions_(other.roptions_),
    prefix_(other.prefix_),
    start_time_(other.start_time_),
    time_limit_(other.time_limit_){
      options_ = new Ipopt::OptionsList();
      *options_ = *other.options_;
  }
  
  TNLPSolver::~TNLPSolver()
  {}


  bool
  TNLPSolver::zeroDimension(const Ipopt::SmartPtr<Ipopt::TNLP>& tnlp, ReturnStatus &optimizationStatus)
  {

    int n,m,dum1, dum2;
    Ipopt::TNLP::IndexStyleEnum dum3;
    tnlp->get_nlp_info(n,m,dum1, dum2, dum3);
    double * x_l = new double[n];
    double * x_u = new double[n];
   
    double * g_l = (m>0) ? new double[m] : NULL;
    double * g_u = (m >0) ? new double[m] : NULL;
    
    tnlp->get_bounds_info(n, x_l, x_u, m, g_l , g_u);
    

    for(int i = 0 ; i < n ; i++) {
      if(x_u[i] - x_l[i] > 1e-5)
	{
	  delete [] x_l;
	  delete [] x_u;
          if(m > 0){
	    delete [] g_l;
	    delete [] g_u;
          }
	  return 0;
	}
    }

    //Problem has no variables just check if the unique solution given by the bounds is
    // feasible or not.
    double obj_value;

    tnlp->eval_f(n, x_l, true, obj_value);

    double * x_sol = new double[n];

      
    IpBlasDcopy(n, x_l, 1, x_sol, 1);

    delete [] x_l;
    delete [] x_u;

    double * g_sol = (m > 0) ? new double [m] : NULL;

    tnlp->eval_g(n, x_sol, true, m, g_sol);
    
    optimizationStatus = solvedOptimal;
    for(int i = 0 ; i < m ; i++) {
      if(g_sol[i] - g_l[i] <  - 1e-07 || g_sol[i] - g_u[i] > 1e-07) {
        optimizationStatus = provenInfeasible;
	
	delete [] g_l;
	delete [] g_u;
	double * lam = (m > 0) ? new double[m]: NULL;
	CoinFillN(lam,m,0.);
	double * z = new double[n];
	CoinFillN(z,n,0.);
	tnlp->finalize_solution(Ipopt::LOCAL_INFEASIBILITY,
			       n, x_sol, NULL, NULL, 
			       m, g_sol, NULL, obj_value, NULL, NULL);
	if (m > 0) delete [] lam;
	delete [] z;
	if (m > 0) delete [] g_sol;
	delete [] x_sol;

        return 1;
      }
    }
    if (m > 0) delete [] g_l;
    if (m > 0) delete [] g_u;

    double * lam = (m > 0) ? new double[m] : NULL;
    CoinFillN(lam,m,0.);
    double * z = new double[n];
    CoinFillN(z,n,0.);
    tnlp->finalize_solution(Ipopt::SUCCESS,
			   n, x_sol, z, z,
			   m, g_sol, lam, obj_value, NULL, NULL);
    if (m > 0) delete [] lam;
    delete [] z;
    if (m > 0) delete [] g_sol;
    delete [] x_sol;
    return 1;
  }

void
TNLPSolver::UnsolvedError::printError(std::ostream &os)
{
  os<<solverName()<<" exited with error code "<<errorNum_<<" "<<errorName()<<std::endl;
}

void
TNLPSolver::UnsolvedError::writeDiffFiles(const std::string prefix) const{
  const int numcols = model_->num_variables();
  const int numrows = model_->num_constraints();
  
  const double * currentLower = model_->x_l();
  const double * currentUpper = model_->x_u();

  const double * originalLower = model_->orig_x_l();
  const double * originalUpper = model_->orig_x_u();
  CoinRelFltEq eq;
  std::string fBoundsName = prefix + name_;
  fBoundsName+="_bounds";
  
  std::string fModName = fBoundsName + ".mod";
  std::ofstream fBounds;
  std::ofstream fMod;

  /** Reader variables names.*/
  bool hasVarNames = 0;
  NamesReader reader(name_,".col");
  
  if(reader.readFile())
      hasVarNames=1;
  if(hasVarNames)
    fMod.open(fModName.c_str());
  fBounds.open(fBoundsName.c_str());
    
  for(int i = 0 ; i < numcols ; i++)
    {    
    if(!eq(currentLower[i],originalLower[i]))
      {
        if(hasVarNames)
          fMod<<"bounds"<<i<<": "
	      <<reader.name(i)<<" >= "
	      <<currentLower[i]<<";\n";


	fBounds<<"LO"<<"\t"<<i<<"\t"<<currentLower[i]<<std::endl;
    }
    if(!eq(currentUpper[i],originalUpper[i]))
      {
	if(hasVarNames)
	  fMod<<"bounds"<<i<<": "
	      <<reader.name(i)<<" <= "
	      <<currentUpper[i]<<";\n";
	
        fBounds<<"UP"<<"\t"<<i<<"\t"<<currentUpper[i]<<std::endl;
      }
    }
  
    //write a file with starting point
    std::string fStartPointName = name_;
    fStartPointName+="_start";



    const double * primals = model_->x_init();
    const double * duals = model_->duals_init();

    if(!primals)//No starting point no output
      {
	std::cerr<<"A failure has occured but no starting point exists"<<std::endl;
	return;
      }

    std::ofstream fStartPoint(fStartPointName.c_str());
    fStartPoint.precision(17);
    fStartPoint<<numcols<<"\t"<<2*numcols+numrows<<std::endl;
    for(int i = 0 ; i < numcols ; i++)
    fStartPoint<<primals[i]<<std::endl;
    int end = 2*numcols + numrows;
    if(duals)
      {
	for(int i = 0 ; i < end; i++)
	  fStartPoint<<duals[i]<<std::endl;
      }

 
}

  /** Say if an optimization status for a problem which failed is recoverable
      (problem may be solvable).*/
  bool 
  TNLPSolver::isRecoverable(ReturnStatus &r){
    return (r >=0 || (r != illDefinedProblem && r != notEnoughFreedom && r != illegalOption && r != computationError && r != timeLimit) );
  }

/** Initialize the options and the journalist.*/
void 
TNLPSolver::initializeOptionsAndJournalist(){
  prefix_ = "bonmin.";
  options_ = new Ipopt::OptionsList();
  
  journalist_= new Ipopt::Journalist();
  roptions_ = new Bonmin::RegisteredOptions();
  
  try{
    Ipopt::SmartPtr<Ipopt::Journal> stdout_journal =
    journalist_->AddFileJournal("console", "stdout", Ipopt::J_ITERSUMMARY);
    
    options_->SetJournalist(journalist_);
    options_->SetRegisteredOptions(GetRawPtr(roptions_));
  }
  catch (Ipopt::IpoptException &E){
    E.ReportException(*journalist_);
    throw E;
  }
  catch(std::bad_alloc){
    journalist_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN, "\n Not enough memory .... EXIT\n");
    throw CoinError("TNLPSolver", "initializeOptionsAndJournalist", "Not enough memory");
  }
#ifndef NO_CATCH_ALL
  catch(...){
    Ipopt::IpoptException E("Uncaught exception in FilterSolver::FilterSolver()",
                            "BonFilterSolver.cpp",-1);
    throw E;
  }
#endif
  
}


}//end namespace Bonmin


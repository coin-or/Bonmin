#include "BonTNLPSolver.hpp"
#include "IpBlas.hpp"
namespace Bonmin{
  using namespace Ipopt;

  TNLPSolver::TNLPSolver()
  {
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
    double * g_l = new double[m];
    double * g_u = new double[m];
    
    tnlp->get_bounds_info(n, x_l, x_u, m, g_l , g_u);
    

    for(int i = 0 ; i < n ; i++) {
      if(x_u[i] - x_l[i] > 1e-5)
	{
	  delete [] x_l;
	  delete [] x_u;
	  delete [] g_l;
	  delete [] g_u;
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

    double * g_sol = new double [m];

    tnlp->eval_g(n, x_sol, true, m, g_sol);
    
    optimizationStatus = solvedOptimal;
    for(int i = 0 ; i < m ; i++) {
      if(g_sol[i] - g_l[i] <  - 1e-07 || g_sol[i] - g_u[i] > 1e-07) {
        optimizationStatus = provenInfeasible;
	
	delete [] g_l;
	delete [] g_u;
	tnlp->finalize_solution(Ipopt::LOCAL_INFEASIBILITY,
			       n, x_sol, NULL, NULL, 
			       m, g_sol, NULL, obj_value);
	delete [] g_sol;
	delete [] x_sol;

        return 1;
      }
    }
    delete [] g_l;
    delete [] g_u;
    tnlp->finalize_solution(Ipopt::SUCCESS,
			   n, x_sol, NULL, NULL,
			   m, g_sol, NULL, obj_value);
    delete [] g_sol;
    delete [] x_sol;
    return 1;
  }

void
TNLPSolver::UnsolvedError::printError(std::ostream &os)
{
  os<<solverName()<<" exited with error code "<<errorNum_<<" "<<errorName()<<std::endl;
}


}//end namespace Bonmin


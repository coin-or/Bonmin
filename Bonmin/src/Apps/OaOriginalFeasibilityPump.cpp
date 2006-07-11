// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 06/18/2005

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <iomanip>


// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiCbcSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcBranchUser.hpp"
#include "CbcCompareUser.hpp"
#include "CbcCompareActual.hpp"
#include "CbcCutGenerator.hpp"
//#include "CbcHeuristicUser.hpp"
#include "SimpleIpoptInterface.hpp"
#include "IpCbcDummyHeuristic.hpp"
#include "IpCbcOACutGenerator.hpp"
#include "IpCbcOACutGenerator2.hpp"

#include "AmplTMINLP.hpp"

#include "CglGomory.hpp"
//#include "CglProbing.hpp"

//#include "ClpQuadInterface.hpp"

// Heuristics would need adapting

//#include "CbcHeuristic.hpp"


// Time
#include "CoinTime.hpp"


int main (int argc, char *argv[])
{

  // Define a Solver which inherits from OsiClpsolverInterface -> OsiSolverInterface

  using namespace Ipopt;

  SmartPtr<IpoptApplication> app = new IpoptApplication();

  app->Jnlst()->Printf(J_ERROR, J_MAIN, "\n\n\n*************************************************************\n");
  app->Jnlst()->Printf(J_ERROR, J_MAIN, "*** Running minlp with AMPL Model  **************************\n");
  app->Jnlst()->Printf(J_ERROR, J_MAIN, "*************************************************************\n\n\n");

  // Read in model using argv[1]
  char * pbName = new char[strlen(argv[1])+1];
  strcpy(pbName, argv[1]);
  SmartPtr<TMINLP> ampl_tminlp = new AmplTMINLP(app->Jnlst(), argv);
  SimpleIpoptInterface solver1(app,ampl_tminlp);

  solver1.initialSolve();
  int nMajorIt = 0;
  bool solved = 0;
  double time=-CoinCpuTime();
  while(!solved && nMajorIt < 50) {
    nMajorIt++;
    const double * colsol = solver1.getColSolution();
    int * inds = new int[solver1.getNumCols()];
    double * x = new double[solver1.getNumCols()];
    int k = 0;
    for(int i = 0 ; i < solver1.getNumCols() ; i++) {
      if(solver1.isInteger(i)) {
        inds[k] = i;
        x[k++] = floor(colsol[i] + 0.5);
        std::cout<<"Var "<<i<<" value "<<x[k-1]<<std::endl;
      }
    }
    OsiCuts cs;
    double dist = solver1.getFeasibilityOuterApproximation( k, x, inds, cs);
    if(dist < 1e-06)
      solved = 1;
    //Change the objective of the MIP to get the closest point to rounding of NLP optimum
  }
  time+=CoinCpuTime();
  if(solved)
    std::cout<<pbName<<" Feasible solution found in "<<time<<" seconds, "<<nMajorIt<<" major iterations"<<std::endl;
  else {
    std::cout<<"Problem aborted on iteration limit elapsed time : "<<time<<",  "<<nMajorIt<<" major iterations"<<std::endl;

  }
  delete [] pbName;
  return 0;
}

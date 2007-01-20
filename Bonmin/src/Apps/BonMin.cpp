// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006 
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Andreas Waechter, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 02/15/2006


#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <iomanip>
#include <fstream>

#include "CoinTime.hpp"
#include "BonAmplInterface.hpp"
#include "BonIpoptSolver.hpp"
#include "BonCbc.hpp"

using namespace Bonmin;

int main (int argc, char *argv[])
{
  using namespace Ipopt;
  
  AmplInterface * nlp_and_solver; 
  char * pbName = NULL;

  if(argc > 1) {
    pbName = new char[strlen(argv[1])+1];
    strcpy(pbName, argv[1]);
  }
  else { //will just output usage
    Ipopt::SmartPtr<IpoptSolver> dummy_ipopt = new IpoptSolver;      
    nlp_and_solver = new AmplInterface(argv, GetRawPtr(dummy_ipopt));
    delete nlp_and_solver;
    return 0;
  }

  double time1 = CoinCpuTime();
  try {
    nlp_and_solver = new AmplInterface(argv);
    BonminCbcParam par;
    Bab bb;
    par(nlp_and_solver);
    bb(nlp_and_solver, par);//do branch and bound

    std::cout.precision(10);

    std::cout<<pbName<<" \t";
    std::string message;
    std::string status;
    if(bb.mipStatus()==Bab::FeasibleOptimal) {
      status = "\t\"Finished\"";
      message = "\nbonmin: Optimal";
    }
    else if(bb.mipStatus()==Bab::ProvenInfeasible) {
      status = "\t\"Finished\"";
      message = "\nbonmin: Infeasible problem";
    }
    else if(bb.mipStatus()==Bab::Feasible) {
      status = "\t\"Not finished\"";
      message = "\n Optimization not finished.";
    }
    else if(bb.mipStatus()==Bab::NoSolutionKnown) {
      status = "\t\"Not finished\"";
      message = "\n Optimization not finished.";
    }

  if(0)// To output a line for building tables
    std::cout<<status<<"\t"<<CoinCpuTime()-time1<<"\t"
	     <<bb.bestObj()<<"\t"
	     <<bb.numNodes()<<"\t"
	     <<bb.iterationCount()<<"\t"
	     <<nlp_and_solver->totalNlpSolveTime()<<"\t"
	     <<nlp_and_solver->nCallOptimizeTNLP()<<"\t"
	     <<std::endl;
    nlp_and_solver->writeAmplSolFile(message,bb.bestSolution(),NULL);

  }
  catch(TNLPSolver::UnsolvedError *E) {
     E->writeDiffFiles();
     E->printError(std::cerr);
    //There has been a failure to solve a problem with Ipopt.
    //And we will output file with information on what has been changed in the problem to make it fail.
    //Now depending on what algorithm has been called (B-BB or other) the failed problem may be at different place.
    //    const OsiSolverInterface &si1 = (algo > 0) ? nlpSolver : *model.solver();
  }
  catch(OsiTMINLPInterface::SimpleError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
	     <<std::endl
	     <<E.message()<<std::endl;
  }
  catch(CoinError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
	     <<std::endl
	     <<E.message()<<std::endl;
  }
  catch (Ipopt::OPTION_INVALID &E)
  {
   std::cerr<<"Ipopt exception : "<<E.Message()<<std::endl;
  }
  catch(...) {
    std::cerr<<pbName<<" unrecognized exception"<<std::endl;
    std::cerr<<pbName<<"\t Finished \t exception"<<std::endl;
    throw;
  }

  delete [] pbName;
  delete nlp_and_solver;
  return 0;
}


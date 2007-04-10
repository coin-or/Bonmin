// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pietro Belloti, Carnegie Mellon University,
// Pierre Bonami, International Business Machines Corporation
//
// Date : 12/19/2006


#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <iomanip>
#include <fstream>

#include "CoinTime.hpp"
#include "BonminConfig.h"
#include "BonCouenneInterface.hpp"
#include "BonIpoptSolver.hpp"
#include "BonCouenneCbc.hpp"

#ifdef COIN_HAS_FILTERSQP
#include "BonFilterSolver.hpp"
#endif
namespace Bonmin{
extern int usingCouenne;}
using namespace Bonmin;

int main (int argc, char *argv[])
{
  using namespace Ipopt;
  Bonmin::usingCouenne = 1;  
  CouenneInterface * nlp_and_solver; 

  // We need to build dummy solver objects to get the options,
  // determine which is the solver to use and register all the options
  Ipopt::SmartPtr<IpoptSolver> dummy_ipopt = new IpoptSolver;
  OsiTMINLPInterface forOption(GetRawPtr(dummy_ipopt));

  int solverUsed = 0; // 0 is Ipopt, 1 is Filter
  forOption.solver()->Options()->GetEnumValue("nlp_solver", solverUsed,"bonmin.");

  char * pbName = NULL;
  if(argc > 1)
  {
    pbName = new char[strlen(argv[1])+1];
    strcpy(pbName, argv[1]);
  }
  else //will just output usage
  {
    Ipopt::SmartPtr<IpoptSolver> ipoptSolver = new IpoptSolver;
    nlp_and_solver = new CouenneInterface (argv, GetRawPtr (ipoptSolver));
    delete nlp_and_solver;
    return 0;
  }
  double time1 = CoinCpuTime();
  try {
  Ipopt::SmartPtr<TNLPSolver> solver;
  if(solverUsed == 0)
    solver = new IpoptSolver;
  else if(solverUsed == 1)
#ifdef COIN_HAS_FILTERSQP
    solver = new FilterSolver;
#else
    {
      std::cerr<<"filterSQP is not properly configured for using into Bonmin"<<std::endl
               <<"be sure to run the configure script with options:"<<std::endl
               <<"--with-filtersqp_lib=\"<path_to_filter_library>\""<<std::endl
               <<"--with-filtersqp_incdir=\"\""<<std::endl;
               throw -1;
      }
#endif
  else
    {
      std::cerr<<"Trying to use unknown solver."<<std::endl;
    }
  nlp_and_solver = new CouenneInterface(argv, solver);
    BonminCbcParam par;
    CouenneBab bb;
    // Eventually change some default in a custom application
    nlp_and_solver->setAppDefaultOptions(nlp_and_solver->solver()->Options());

    par(nlp_and_solver);
    bb(nlp_and_solver, par);//do branch and bound

    std::cout.precision(10);

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

    if (0) {// print statistics in LaTeX format

      char *basename = strrchr (pbName, '/');
      if (!basename) basename = pbName;
      else basename++;

      printf (" %-25s & %8.2f &", basename, CoinCpuTime () - time1);

      if (fabs (bb.bestBound()) < 1e40) 
	printf (" %12.3f &", bb.bestBound());
      else printf (" %8s     &", "inf_dual");
	  
      if (fabs (bb.bestObj()) < 1e40) 
	printf (" %12.3f &", bb.bestObj());
      else printf (" %8s     &", "inf_prim");
	  
      printf ("%7d & %7d & %-20s\\\\\n ",
	      bb.numNodes(),
	      bb.iterationCount(),
	      //	      nlp_and_solver->totalNlpSolveTime(),
	      //	      nlp_and_solver->nCallOptimizeTNLP(),
	      status.c_str());
      
    }

    if (0) { // print statistics
      std::cout<<pbName << status
	       <<CoinCpuTime()-time1<<"\t"
	       <<bb.bestBound()<<"\t"
	       <<bb.bestObj()<<"\t"
	       <<bb.numNodes()<<"\t"
	       <<bb.iterationCount()<<"\t"
	       <<nlp_and_solver->totalNlpSolveTime()<<"\t"
	       <<nlp_and_solver->nCallOptimizeTNLP()<<"\t"
	       <<status
	       <<std::endl;
    }

    //  nlp_and_solver->writeAmplSolFile(message,bb.bestSolution(),NULL);

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
//  catch(...) {
//    std::cerr<<pbName<<" unrecognized excpetion"<<std::endl;
//    std::cerr<<pbName<<"\t Finished \t exception"<<std::endl;
//    throw;
//  }

  delete [] pbName;
  delete nlp_and_solver;
  return 0;
}


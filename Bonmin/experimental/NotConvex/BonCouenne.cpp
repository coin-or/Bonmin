// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pietro Belotti, Carnegie Mellon University,
// Pierre Bonami, International Business Machines Corporation
//
// Date : 12/19/2006


#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <iomanip>
#include <fstream>

#include <stdlib.h>

#include "CoinTime.hpp"
#include "BonminConfig.h"
#include "BonCouenneInterface.hpp"
#include "BonIpoptSolver.hpp"

#include "BonCouenneSetup.hpp"

#include "BonCbc.hpp"
#ifdef COIN_HAS_FILTERSQP
#include "BonFilterSolver.hpp"
#endif

#include "CbcCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "CouenneCutGenerator.hpp"

using namespace Bonmin;

// global variable (for testing purposes only!) /////////////
//double fakeCutOff = +1e50;

int main (int argc, char *argv[])
{
  using namespace Ipopt;

  char * pbName = NULL;
  if(argc > 1)
  {
    pbName = new char[strlen(argv[1])+1];
    strcpy(pbName, argv[1]);
  }

  /*if (argc > 2) {
    fakeCutOff = atof (argv [2]);
    printf ("cutoff = %g\n", fakeCutOff);
    argc = 2;
    *(argv [2]) = 0;
    }*/

  double time1 = CoinCpuTime();

  try {
    CouenneSetup bonmin;
    bonmin.InitializeCouenne (argv);
    Bab bb;

#if 0
    CouenneFeasibility feasibility;
    bb.model().setProblemFeasibility (feasibility);
#endif

    //////////////////////////////////

    bb (bonmin); // do branch and bound

    //    double opt_x0 = bb.model () . bestSolution () [0];

    /*    printf ("#m=1,S=0 # brtest\n\
%10g %10g # brtest\n\
%10g %10g # brtest\n\
 # brtest\n\
#m=-1,S=0 # brtest\n\
# brtest\n", opt_x0, log (opt_x0) - 1, opt_x0, log (opt_x0) + 1);*/

    //////////////////////////////////

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

    if (0) { // print statistics in LaTeX format

      ////////////////////////////////
      int nr=-1, nt=-1;
      double st=-1;

      // CAUTION: assuming first cut generator is our CouenneCutGenerator

      CouenneCutGenerator *cg = NULL;

      if (bb.model (). cutGenerators ())
	cg = dynamic_cast <CouenneCutGenerator *> 
	  (bb.model (). cutGenerators () [0] -> generator ());

      if (cg) cg -> getStats (nr, nt, st);
      else printf ("Warning, could not get pointer to CouenneCutGenerator\n");

      char *basename = strrchr (pbName, '/');
      if (!basename) basename = pbName;
      else basename++;

      CouenneProblem *cp = cg ? cg -> Problem () : NULL;

      if (cg && !cp) printf ("Warning, could not get pointer to problem\n");

      printf ("::: %-15s & %6d & %6d & %6d & %6d & %10d & %10d & %8.3f & ", 
	      basename,
	      (cp) ? cp -> nOrig     () : -1, 
	      (cp) ? cp -> nIntVars  () : -1, 
	      (cp) ? cp -> nOrigCons () : -1,
	      (cp) ? (cp -> nVars   () - 
		      cp -> nOrig   ()): -1,
	      nr, nt, st);

      /////////////////////////////////

      double timeLimit = 0, obj = bb.model (). getObjValue ();
      bonmin.options() -> GetNumericValue ("time_limit", timeLimit, "bonmin.");

      if (CoinCpuTime () - time1 > timeLimit) {

	// time limit reached, print upper and (in brackets) lower

	double obj = bb.model (). getObjValue ();

	if (fabs (obj) < 9e12) 
	  printf    (" %18.9f &", bb.bestObj ());
	else printf (" %8s     &", "inf_prim");

	if (fabs (bb.bestBound()) < 9e12) 
	  printf    (" (%18.9f) &", bb.bestBound ());
	else printf (" %8s       &", "inf_dual");
      }
      else {
	// time limit not reached, print upper and time

	if (fabs (obj) < 9e12) 
	  printf    (" %18.9f &", bb.bestObj ());
	else printf (" %8s     &", "inf_prim");
	  
	if (fabs (bb.bestBound()) < 9e12) 
	  printf    (" %12.3f   &", CoinCpuTime () - time1);
	else printf (" %8s       &", "inf_dual");
      }

      printf ("%7d & %7d \\\\\n",
	      bb.numNodes(),
	      bb.iterationCount());
	      //	      nlp_and_solver->totalNlpSolveTime(),
	      //	      nlp_and_solver->nCallOptimizeTNLP(),
      //	      status.c_str());
    }

//    nlp_and_solver -> writeAmplSolFile (message, bb.bestSolution (), NULL);
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
  return 0;
}

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

int main (int argc, char *argv[])
{
  using namespace Ipopt;

  char * pbName = NULL;
  double time1 = CoinCpuTime();

  try {

    Bab bb;
    bb.setUsingCouenne (true);

    CouenneSetup bonmin;
    bonmin.InitializeCouenne (argv);

#if 0
    CouenneFeasibility feasibility;
    bb.model().setProblemFeasibility (feasibility);
#endif

    //////////////////////////////////

    bb (bonmin); // do branch and bound

    //////////////////////////////////

    std::cout.precision (10);

    /*std::string message, status;

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
      }*/

    if (bonmin.displayStats ()) { // print statistics

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

      CouenneProblem *cp = cg ? cg -> Problem () : NULL;

      if (cg && !cp) printf ("Warning, could not get pointer to problem\n");
      else
	printf ("Stats: %-15s %4d [var] %4d [int] %4d [con] %4d [aux] "
		"%6d [root] %8d [tot] %6g [sep] %8g [time] %8g [bb] "
		"%20e [lower] %20e [upper] %7d [nodes]\n",// %s %s\n",
		cp ? cp -> problemName ().c_str () : "unknown",
		(cp) ? cp -> nOrig     () : -1, 
		(cp) ? cp -> nIntVars  () : -1, 
		(cp) ? cp -> nOrigCons () : -1,
		(cp) ? (cp -> nVars   () - 
			cp -> nOrig   ()): -1,
		nr, nt, st, 
		CoinCpuTime () - time1,
		cg ? (CoinCpuTime () - cg -> rootTime ()) : - CoinCpuTime (),
		bb.bestBound (),
		//bb.bestObj (),
		bb.model (). getObjValue (),
		bb.numNodes ());
		//bb.iterationCount ());
		//status.c_str (), message.c_str ());

      /*
      /////////////////////////////////

      double timeLimit = 0, obj = bb.model (). getObjValue ();
      bonmin.options() -> GetNumericValue ("time_limit", timeLimit, "bonmin.");

      if (CoinCpuTime () - time1 > timeLimit) {

	// time limit reached, print upper and (in brackets) lower

	double obj = bb.model (). getObjValue ();

	if (fabs (obj) < 9e12) 
	  printf    ("%18.9f &", bb.bestObj ());
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
      */
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

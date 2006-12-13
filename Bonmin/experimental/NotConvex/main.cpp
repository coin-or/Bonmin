// (C) Copyright International Business Machines Corporation (IBM and Carnegie Mellon University 2006 
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 02/15/2006


#include "CouenneProblem.h"
#include <iomanip>
#include <fstream>
#include <iostream>

#include "CoinTime.hpp"
#include "BonminConfig.h"
#include "BonAmplInterface.hpp"
#include "BonIpoptSolver.hpp"
#include "BonCbc.hpp"

#ifdef COIN_HAS_FILTERSQP
#include "BonFilterSolver.hpp"
#endif

#include <unistd.h>

#include <stdio.h>
#include <iostream>

#include <math.h>

/* AMPL includes */
#include "asl.h"
#include "asl_pfgh.h"
#include "getstub.h"


int main (int argc, char **argv) {

// Some Bonmin stuff

  //Need to build a dummy solver to read options
  Ipopt::SmartPtr<Bonmin::IpoptSolver> dummy_ipopt = 
                               new Bonmin::IpoptSolver;
  dummy_ipopt->Initialize("bonmin.opt");
  int solverUsed = 0; //O if it is ipopt 1 if it is filterSqp

  dummy_ipopt->Options()->GetEnumValue("nlp_solver", solverUsed , "bonmin.");

    char * pbName = NULL;
  if(argc > 1)
  {
    pbName = new char[strlen(argv[1])+1];
    strcpy(pbName, argv[1]);
  }
  else //will just output usage
  {
    Ipopt::SmartPtr<Bonmin::IpoptSolver> ipoptSolver = new Bonmin::IpoptSolver;
    Bonmin::AmplInterface * nlp_and_solver = new Bonmin::AmplInterface(argv, GetRawPtr(ipoptSolver));
    delete nlp_and_solver;
    return 0;
  }

  try{
      Ipopt::SmartPtr<Bonmin::TNLPSolver> solver;
  if(solverUsed == 0)
    solver = new Bonmin::IpoptSolver;
  else if(solverUsed == 1)
#ifdef COIN_HAS_FSQP
    solver = new Bonmin::FilterSolver;
#else
    {
      std::cerr<<"filterSQP is not propoerly configured for using into Bonmin"<<std::endl
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
   Bonmin::AmplInterface * nlp_and_solver = new Bonmin::AmplInterface(argv, solver);


//Now we come to the real meat

  // Declare some Couenne problem. 
  CouenneProblem p;
  //Get the ASL pointer and pass it to Pietro
  const ASL_pfgh* asl = nlp_and_solver->amplModel()->AmplSolverObject();
//  p.readNl(ASLPointer);



    delete nlp_and_solver;
  }
  catch(Bonmin::TNLPSolver::UnsolvedError *E)
  {
    E->writeDiffFiles();
    E->printError(std::cerr);
  }
  catch(Bonmin::OsiTMINLPInterface::SimpleError &E) {
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
    std::cerr<<pbName<<" unrecognized excpetion"<<std::endl;
    std::cerr<<pbName<<"\t Finished \t exception"<<std::endl;
    throw;
  }

  delete [] pbName;
  return 0;


} 

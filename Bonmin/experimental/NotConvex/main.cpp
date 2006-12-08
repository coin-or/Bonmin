// (C) Copyright International Business Machines Corporation (IBM and Carnegie Mellon University 2006 
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Andreas Waechter, International Business Machines Corporation
// Pierre Bonami, International Business Machines Corporation
// Pietro Belotti, Carnegie Mellon University
//
// Date : 02/15/2006


#include "problem.h"
#include <iomanip>
#include <fstream>
#include <iostream>

#include "CoinTime.hpp"
#include "BonminConfig.h"
#include "BonAmplInterface.hpp"
#include "BonIpoptSolver.hpp"
#include "BonCbc.hpp"

#ifdef COIN_HAS_FSQP
#include "BonFilterSolver.hpp"
#endif

#include <unistd.h>

#include <stdio.h>
#include <iostream>

#include <math.h>



int main (int argc, char **argv) {

// Some Bonmin stuff

  //Need to build a dummy solver to read options
  Ipopt::SmartPtr<Bonmin::IpoptSolver> dummy_ipopt = 
                               new Bonmin::IpoptSolver;
  dummy_ipopt->Initialize("bonmin.opt");
  int solverUsed = 0; //O if it is ipopt 1 if it is filterSqp

  dummy_ipopt()->Options()->GetEnumValue("nlp_solver", solverUsed , "bonmin.");

    char * pbName = NULL;
  if(argc > 1)
  {
    pbName = new char[strlen(argv[1])+1];
    strcpy(pbName, argv[1]);
  }
  else //will just output usage
  {
    Ipopt::SmartPtr<IpoptSolver> ipoptSolver = new IpoptSolver;
    nlp_and_solver = new AmplInterface(argv, GetRawPtr(ipoptSolver));
    delete nlp_and_solver;
    return 0;
  }

  try{
      Ipopt::SmartPtr<TNLPSolver> solver;
  if(solverUsed == 0)
    solver = new IpoptSolver;
  else if(solverUsed == 1)
#ifdef COIN_HAS_FSQP
    solver = new FilterSolver;
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
   Bonmin::AmplInterface * nlp_and_solver = new AmplInterface(argv, solver);


//Now we come to the real meat
  Problem p;
  //Get the ASL pointer and pass it to Pietro
  ASL_pfgh* nlp_and_solver->amplModel()->AmplSolverObject();
//  p.init(ASLPointer);




  }
  catch(Bonmin::UnsolvedError *E)
  {
    E->writeDiffFiles();
    E->printError(std::cerr);
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
    std::cerr<<pbName<<" unrecognized excpetion"<<std::endl;
    std::cerr<<pbName<<"\t Finished \t exception"<<std::endl;
    throw;
  }

  delete [] pbName;
  delete nlp_and_solver;
  return 0;


} 

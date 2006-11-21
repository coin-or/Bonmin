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

#include <cassert>
#include <iomanip>

#include "BonminConfig.h"

#include "BonAmplInterface.hpp"
#include "BonIpoptSolver.hpp"
#include "BonBoundsReader.hpp"
#include "BonStartPointReader.hpp"
#ifdef COIN_HAS_FSQP
#include "BonFilterSolver.cpp"
#endif
#include "CoinTime.hpp"

/************************************************************************
 
This mains is used for resolving the problem with fixed bounds and eventually a starting point
 
************************************************************************/

int main (int argc, char *argv[])
{

  using namespace Ipopt;
  using namespace Bonmin;

  // Read in model using argv[1]
  char * pbName = new char[strlen(argv[1])+1];
  strcpy(pbName, argv[1]);
  std::string nodeFileName;
  if(argc > 2)
    nodeFileName=argv[2];
  std::string startingPointFile ="";
  if(argc>3)
    startingPointFile = argv[3];

  //Give ampl an argv which doesn't crash him.
  char ** myArgv = new char *[3];
  myArgv[0]=new char[strlen(argv[0])+1];
  strcpy(myArgv[0],argv[0]);
  myArgv[1]=new char[strlen(argv[1])+1];
  strcpy(myArgv[1],argv[1]);
  myArgv[2]= NULL;//new char[1];

  //We need to build dummy solver objects to get the options, determine which is the solver to use and register all the options
  Ipopt::SmartPtr<IpoptSolver> dummy_ipopt = new IpoptSolver;
  OsiTMINLPInterface forOption(GetRawPtr(dummy_ipopt));


  int solverUsed = 0; // 0 is Ipopt, 1 is Filter
  forOption.solver()->Options()->GetEnumValue("nlp_solver", solverUsed,"bonmin.");


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

  Bonmin::AmplInterface nlpSolver(myArgv, solver);

  Ipopt::SmartPtr<Ipopt::OptionsList> Options =
    nlpSolver.retrieve_options();

  nlpSolver.messageHandler()->setLogLevel(2);

  try
    {
      std::cout<<nodeFileName<<std::endl;
      // Read the bounds and change them in Ipopt
      if(argc>2) {
	Bonmin::BoundsReader bounds(nodeFileName);
	bounds.readAndApply(&nlpSolver);
      }
      if(argc>3) {
	Bonmin::StartPointReader init(startingPointFile);
	init.readAndApply(&nlpSolver);
      }
      
      nlpSolver.solver()->turnOnOutput();
      nlpSolver.initialSolve();
      
      //Print out integer variable values
      for(int i = 0 ; i <nlpSolver.getNumCols() ; i++) {
	if (nlpSolver.isInteger(i)) {
	  std::cout<<"x[ "<<i<<"] = "<<nlpSolver.getColSolution()[i]<<std::endl;
	}
      }
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
  catch(...) {
    std::cerr<<pbName<<" unrecognized excpetion"<<std::endl;
    std::cerr<<pbName<<"\t Finished \t exception"<<std::endl;
    throw;
  }
  
  
  delete [] pbName;
  delete [] myArgv[0];
  delete [] myArgv[1];
  delete [] myArgv;
  return 0;
}

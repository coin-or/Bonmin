// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006 
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Andreas Waechter, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 02/15/2006


#include <cassert>
#include <iomanip>

#include "BonminConfig.h"

#include "CoinPragma.hpp"
#include "BonAmplInterface.hpp"
#include "BonIpoptSolver.hpp"
#include "BonBoundsReader.hpp"
#include "BonStartPointReader.hpp"
#include "CoinTime.hpp"

#include "BonAmplSetup.hpp"


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


    BonminAmplSetup bonmin;
    bonmin.initialize(myArgv);
    Bonmin::OsiTMINLPInterface& nlpSolver = *bonmin.nonlinearSolver();
  
    Ipopt::SmartPtr<Ipopt::OptionsList> Options =
      nlpSolver.options();

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
      
      nlpSolver.solver()->forceSolverOutput(4);
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

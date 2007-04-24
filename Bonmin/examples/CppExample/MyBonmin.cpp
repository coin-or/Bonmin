// (C) Copyright Carnegie Mellon University 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  03/17/2006


#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <iomanip>
#include <fstream>

#include "CoinTime.hpp"

#include "BonOsiTMINLPInterface.hpp"
#include "BonIpoptSolver.hpp"
#include "MyTMINLP.hpp"
#include "BonCbc2.hpp"
#include "BonBonminSetup.hpp"

#include "BonOACutGenerator2.hpp"
#include "BonEcpCuts.hpp"
#include "BonOaNlpOptim.hpp"


int main (int argc, char *argv[])
{
  using namespace Ipopt;
  using namespace Bonmin;
  SmartPtr<MyTMINLP> tminlp = new MyTMINLP;
  
  BasicSetup b;
  
  //Register an additional option
  b.roptions()->AddStringOption2("print_solution","Do we print the solution or not?",
                                 "yes",
                                 "no", "No, we don't.",
                                 "yes", "Yes, we do.",
                                 "A longer comment can be put here");
  
  
  // Register all the bonmin options.
  BonminSetup::registerAllOptions(b.roptions());
  
  // Here we can change the default value of some Bonmin or Ipopt option
  b.options()->SetNumericValue("bonmin.time_limit", 5); //changes bonmin's time limit
  b.options()->SetStringValue("mu_oracle","loqo");
  
  //Here we can read one or several option files
  b.Initialize("Mybonmin.opt");
  b.Initialize("bonmin.opt");

  
  // Now we can obtain the value of the new option
  int printSolution;
  b.options()->GetEnumValue("print_solution", printSolution,"");
  if(printSolution == 1){
    tminlp->printSolutionAtEndOfAlgorithm();
  }
  
  BonminSetup bonmin;
  bonmin.setBasicOptions(b);
  bonmin.initializeBonmin(GetRawPtr(tminlp));



  //Set up done, now let's branch and bound
  double time1 = CoinCpuTime();
  try {
    Bab2 bb;
    bb(bonmin);//process parameter file using Ipopt and do branch and bound


  }
  catch(TNLPSolver::UnsolvedError *E) {
    //There has been a failure to solve a problem with Ipopt.
    std::cerr<<"Ipopt has failed to solve a problem"<<std::endl;
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


  return 0;
}


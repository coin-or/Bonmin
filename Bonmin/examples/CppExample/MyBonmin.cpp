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
  SmartPtr<TMINLP> tminlp = new MyTMINLP;
  BonminSetup bonmin;
  bonmin.initializeBonmin(tminlp);

  bonmin.options()->SetNumericValue("bonmin.time_limit", 1); //changes bonmin's time limit
  bonmin.options()->SetStringValue("mu_oracle","loqo");

  // we can also try and read an option file (this can eventually change options set before, option file always have priority)
//  bonmin.initialize("My_bonmin.opt");

  //Set up done, now let's branch and bound
  double time1 = CoinCpuTime();
  try {
    Bab2 bb;
    bb(bonmin);//process parameter file using Ipopt and do branch and bound

    std::cout.precision(10);

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


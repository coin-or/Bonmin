// (C) Copyright Carnegie Mellon University 2006
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
#include "BonCbc.hpp"



int main (int argc, char *argv[])
{
  using namespace Ipopt;
  using namespace Bonmin;
  SmartPtr<TMINLP> tminlp = new MyTMINLP;
  OsiTMINLPInterface nlpSolver(tminlp, new IpoptSolver);
  
  //Option can be set here directly either to bonmin or ipopt
  nlpSolver.retrieve_options()->SetNumericValue("bonmin.time_limit", 1); //changes bonmin's time limit
  nlpSolver.retrieve_options()->SetStringValue("mu_oracle","loqo");

  // we can also try and read an option file (can eventually change options set before, option file always have priority)
  nlpSolver.readOptionFile("My_bonmin.opt");

  //Set up done, now let's branch and bound
  double time1 = CoinCpuTime();
  try {
    BonminCbcParam par;
    Bab bb;
    par(&nlpSolver);

    bb(&nlpSolver, par);//process parameter file using Ipopt and do branch and bound

    std::cout.precision(10);

    std::string message;
    if(bb.mipStatus()==Bab::FeasibleOptimal) {
      std::cout<<"\t\"Finished\"\t";
      message = "\nbonmin: Optimal";
    }
    else if(bb.mipStatus()==Bab::ProvenInfeasible) {
      std::cout<<"\t\"Finished\"\t";
      message = "\nbonmin: Infeasible problem";
    }
    else if(bb.mipStatus()==Bab::Feasible) {
      std::cout<<"\t\"Not finished\""<<"\t";
      message = "\n Optimization not finished.";
    }
    else if(bb.mipStatus()==Bab::NoSolutionKnown) {
      std::cout<<"\t\"Not finished\""<<"\t";
      message = "\n Optimization not finished.";
    }
    std::cout<<CoinCpuTime()-time1<<"\t"
    <<bb.bestObj()<<"\t"
    <<bb.numNodes()<<"\t"
    <<bb.iterationCount()<<"\t"
    <<std::endl;

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


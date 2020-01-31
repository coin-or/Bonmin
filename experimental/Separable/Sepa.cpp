// (C) Copyright CNRS and others 2010
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, Université de la Méditérannée
// Hassan Hijazi, Orange Labs
//
// Date : 05/22/2010

#include <iomanip>
#include <fstream>

#include "CoinPragma.hpp"
#include "BonminConfig.h"
#include "IpoptConfig.h"
#include "CbcConfig.h"
#include "CoinTime.hpp"
#include "SepaSetup.hpp"
#include "BonCbc.hpp"

#define CATCH_ERRORS
int main (int argc, char *argv[])
{
  using namespace Ipopt;
  char * pbName = NULL;
  
  std::cout<<"Bonmin "
           <<BONMIN_VERSION; 
  std::cout<<" using Cbc "
         <<CBC_VERSION; 
  std::cout<<" and Ipopt "
         <<IPOPT_VERSION<<std::endl; 
  if(argc > 1) {
    pbName = new char[strlen(argv[1])+1];
    strcpy(pbName, argv[1]);
  }
#define CATCH_ERRORS 
#ifdef CATCH_ERRORS
  try
#endif 
  {

    //FILE * fp = fopen("log","w");
    Sepa::SepaSetup bonmin;
    bonmin.initialize(argv);
    Bonmin::Bab bb;

    bb(bonmin);//do branch and bound

  }
#ifdef CATCH_ERRORS
  catch(Bonmin::TNLPSolver::UnsolvedError *E) {
    E->writeDiffFiles();
    E->printError(std::cerr);
    //There has been a failure to solve a problem with Ipopt.
    //And we will output file with information on what has been changed in the problem to make it fail.
    //Now depending on what algorithm has been called (B-BB or other) the failed problem may be at different place.
    //    const OsiSolverInterface &si1 = (algo > 0) ? nlpSolver : *model.solver();
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
  catch (Ipopt::IpoptException &E)
  {
    std::cerr<<"Ipopt exception : "<<E.Message()<<std::endl;
  }
  catch(...) {
    std::cerr<<pbName<<" unrecognized exception"<<std::endl;
    std::cerr<<pbName<<"\t Finished \t exception"<<std::endl;
    throw;
  }
#endif
  
  delete [] pbName;
  return 0;
}


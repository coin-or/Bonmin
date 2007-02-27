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

#include <iomanip>
#include <fstream>

#include "CoinTime.hpp"

#include "BonminAmplInterface.hpp"
#include "CbcBonmin.hpp"

#ifndef COIN_HAS_ASL 
#error "Ampl Solver Library (ASL) is not properly configured."
#error "Can not compile bonmin executable."
#error "Please check that you have the file asl.h in ThirdPary/ASL subdirectory"
#error "See documentation for instuctions on installing ASL."
#error "And rerun the configure script."
#endif


/** Procedure to ouptut relevant informations in the case of a failure.
    si1 should be the problem solved at a node of the b&b tree, and si2 the original problem.
    Compare the two problems stored in si1 and si2
    and writes files with bounds which have changed.
    Also outputs a file with the starting point of si1.

*/
void writeNodeFiles(const OsiSolverInterface& si1,const IpoptInterface& si2)
{
  const int numcols = si1.getNumCols();
  const int numrows = si1.getNumRows();
  assert( numcols==si2.getNumCols());
  
  const double * currentLower = si1.getColLower();
  const double * currentUpper = si1.getColUpper();

  const double * originalLower = si2.problem()->orig_x_l();
  const double * originalUpper = si2.problem()->orig_x_u();
  CoinRelFltEq eq;
  std::string fBoundsName;
  si2.getStrParam(OsiProbName,fBoundsName);
  fBoundsName+="_bounds";
  
  std::string fModName = fBoundsName + ".mod";
  std::ofstream fBounds;
  std::ofstream fMod;
  bool hasVarNames = 0;
  
  if(si2.getVarNames()!=NULL )
      hasVarNames=1;
  if(hasVarNames)
    fMod.open(fModName.c_str());
  fBounds.open(fBoundsName.c_str());
    
  for(int i = 0 ; i < numcols ; i++)
    {    
    if(!eq(currentLower[i],originalLower[i]))
      {
        if(hasVarNames)
          fMod<<"bounds"<<i<<": "
	      <<si2.getVarNames()[i]<<" >= "
	      <<currentLower[i]<<";\n";


	fBounds<<"LO"<<"\t"<<i<<"\t"<<currentLower[i]<<std::endl;
    }
    if(!eq(currentUpper[i],originalUpper[i]))
      {
	if(hasVarNames)
	  fMod<<"bounds"<<i<<": "
	      <<si2.getVarNames()[i]<<" <= "
	      <<currentUpper[i]<<";\n";
	
        fBounds<<"UP"<<"\t"<<i<<"\t"<<currentUpper[i]<<std::endl;
      }
    }
  
    //write a file with starting point
    std::string fStartPointName;
    si2.getStrParam(OsiProbName,fStartPointName);
    fStartPointName+="_start";



    const IpoptInterface * ipopt = dynamic_cast<const IpoptInterface *>(&si1);
    assert(ipopt);

    const double * primals = ipopt->problem()->x_init();
    const double * duals = ipopt->problem()->duals_init();

    if(!primals)//No starting point no output
      {
	std::cerr<<"A failure has occured but no starting point exists"<<std::endl;
	return;
      }

    std::ofstream fStartPoint(fStartPointName.c_str());
    fStartPoint.precision(17);
    fStartPoint<<numcols<<"\t"<<2*numcols+numrows<<std::endl;
    for(int i = 0 ; i < numcols ; i++)
    fStartPoint<<primals[i]<<std::endl;
    int end = 2*numcols + numrows;
    if(duals)
      {
	for(int i = 0 ; i < end; i++)
	  fStartPoint<<duals[i]<<std::endl;
      }

}

int main (int argc, char *argv[])
{
  using namespace Ipopt;
  
  BonminAmplInterface * nlpSolver; 
  char * pbName = NULL;
  if(argc > 1)
  {
    pbName = new char[strlen(argv[1])+1];
    strcpy(pbName, argv[1]);
  }
  else //will just output usage
  {
    nlpSolver = new BonminAmplInterface(argv);
    delete nlpSolver;
    return 0;
  }
  double time1 = CoinCpuTime();
  try {
  nlpSolver = new BonminAmplInterface(argv, true);
    BonminCbcParam par;
    BonminBB bb;
    par(*nlpSolver);
    //par.fout = fopen("out.txt","w");
    bb(*nlpSolver, par);//do branch and bound

    std::cout.precision(10);
  }
  catch(IpoptInterface::UnsolvedError &E) {
     E.printError(std::cerr);
    //There has been a failure to solve a problem with Ipopt.
    //And we will output file with information on what has been changed in the problem to make it fail.
    //Now depending on what algorithm has been called (B-BB or other) the failed problem may be at different place.
    //    const OsiSolverInterface &si1 = (algo > 0) ? nlpSolver : *model.solver();
#ifdef BONMIN_DEBUG
    writeNodeFiles(*nlpSolver, *nlpSolver);
#endif
  }
  catch(IpoptInterface::SimpleError &E) {
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
  delete nlpSolver;
  return 0;
}


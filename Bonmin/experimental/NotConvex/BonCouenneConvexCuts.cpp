 // (C) Copyright International Business Machines (IBM) 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, International Business Machines
//
// Date :  12/20/2006

#include "BonCouenneConvexCuts.hpp"
#include "BonCouenneInterface.hpp"
#include "CouenneCutGenerator.h"
namespace Bonmin{

double 
CouenneConvCuts::doCouenneConvRounds(OsiSolverInterface &si, 
                     bool leaveSiUnchanged){
  OsiSolverInterface * saveLp = lp_;
  lp_ = &si;
  OsiCuts cs;
  bool saveLeaveSi = leaveSiUnchanged_;
  leaveSiUnchanged_ = leaveSiUnchanged;
  generateCuts(si, cs);
  lp_ = saveLp;
  leaveSiUnchanged_ = saveLeaveSi;
  return objValue_;
}

void
CouenneConvCuts::generateCuts(const OsiSolverInterface &si, 
                      OsiCuts & cs,
                      const CglTreeInfo info) const
{
  double num=CoinDrand48();
  const int & depth = info.level;
  double beta=10000;
  //if(depth == 0) return;
  //   if(num> beta*pow(2,-depth))
  //     return;
  double violation = 1.;//nlp_->getNonLinearitiesViolation(
  //si.getColSolution(), si.getObjValue());
  //  std::cout<<"Constraint violation: "<<violation<<std::endl;

  //Get a random number
  CouenneInterface* ci = dynamic_cast<CouenneInterface *>(nlp_);
  if(violation <= 1e-12)
    return;
  std::cout<<"Generating ECP cuts"<<std::endl;
  solverManip * lpManip = NULL;
  bool infeasible = false;
  const int memory = 5;
  double lastObjs[memory];
  double minImprove = 0.05;
  for(int i = 0 ; i < numRounds_ ; i++)
  {
    const double * colsol = si.getColSolution();
    int numcols = si.getNumCols();
    std::cout<<"Cutting :"<<std::endl;
    for(int kk = 0 ; kk < numcols; kk++){
      std::cout<<"x["<<kk<<"] = "<<colsol[kk]<<"\t";      
    }
    std::cout<<std::endl;
    if( violation > 1e-06)
    {
      int numberCuts =  - cs.sizeRowCuts();
      {

	ci->couenneCg()->updateConv(
				    const_cast<double *>(si.getColSolution()), 
				    const_cast<double *>(si.getColLower()), 
				    const_cast<double *>(si.getColUpper()));
	

	int ncuts = ci->couenneCg()->getncuts();
	for(int i = 0 ; i < ncuts ; i++)
	  {
	    cs.insert(*ci->couenneCg()->getCut(i)); 
	  }
      }
      numberCuts += cs.sizeRowCuts();
      if(numberCuts > 0 && (lp_ || i + 1 < numRounds_ )){
        if(lpManip==NULL) {
          if(lp_ == NULL)
            lpManip = new solverManip(si);
          else
            lpManip = new solverManip(lp_, true,true, 
                                      false,false);
        }
        lpManip->installCuts(cs,numberCuts);
        lpManip->si()->resolve();
	std::cout<<"Cut generation: round "<<i<<", objective "<<lpManip->si()->getObjValue()<<std::endl;
        if(lpManip->si()->isProvenPrimalInfeasible())
        {
          infeasible = true;
	  //      std::cout<<"Stopping CouenneConv generation because problem became infeasible"<<std::endl;
          break;
        }
	if(i >= memory)
	  {
	    double newObj = lpManip->si()->getObjValue();
	    int ind = i % memory;
	    if(fabs((newObj - lastObjs[ind])/newObj)<minImprove)//Is tailing off
	      break;
	    else
	      lastObjs[ind] = newObj;
	  }
	else {
	  lastObjs[i] = lpManip->si()->getObjValue();}
	//        violation =  
	//nlp_->getNonLinearitiesViolation(
	//lpManip->si()->getColSolution(), 
	//lpManip->si()->getObjValue());
	//        std::cout<<"Constraint violation: "<<violation<<std::endl;
         }
      else break;
    }
    else break;
  }
  if(!infeasible){
    if(lpManip != NULL)
      {
	lpManip->si()->resolve();
	if(lpManip->si()->isProvenPrimalInfeasible())
	  objValue_ = 2e50;
	else
	  objValue_ = lpManip->si()->getObjValue();}
  }
  else objValue_ = 2e50;
  if(lpManip)
    {
      if(lp_ != NULL && lpManip != NULL)
	{
	  lpManip->restore();
	}
      
      delete lpManip;
    }
  //  std::cout<<"End ecp cut generation"<<std::endl;
  return;
}

} // end namespace bonmin.

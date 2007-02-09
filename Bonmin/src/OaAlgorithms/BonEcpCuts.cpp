 // (C) Copyright International Business Machines (IBM) 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, International Business Machines
//
// Date :  12/20/2006

#include "BonEcpCuts.hpp"
namespace Bonmin{

double 
EcpCuts::doEcpRounds(OsiSolverInterface &si, 
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
EcpCuts::generateCuts(const OsiSolverInterface &si, 
                      OsiCuts & cs,
                      const CglTreeInfo info) const
{
  double num=CoinDrand48();
  const int & depth = info.level;
  double beta=10000;
  if(depth == 0) return;
  if(num> beta*pow(2,-depth))
    return;
  double violation = nlp_->getNonLinearitiesViolation(
                     si.getColSolution(), si.getObjValue());
  //  std::cout<<"Constraint violation: "<<violation<<std::endl;

  //Get a random number
  if(violation <= 1e-12)
    return;
  std::cout<<"Generating ECP cuts"<<std::endl;
  solverManip * lpManip = NULL;
  bool infeasible = false;
  for(int i = 0 ; i < numRounds_ ; i++)
  {
    if( violation > 1e-06)
    {
      int numberCuts =  - cs.sizeRowCuts();
      const double * toCut = parameter().addOnlyViolated_?
	si.getColSolution():NULL;
      nlp_->getOuterApproximation(cs, 1, toCut, parameter().global_);
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
        if(lpManip->si()->isProvenPrimalInfeasible())
        {
          infeasible = true;
	  //      std::cout<<"Stopping Ecp generation because problem became infeasible"<<std::endl;
          break;
        }
        violation =  nlp_->getNonLinearitiesViolation(
                     lpManip->si()->getColSolution(), 
                     lpManip->si()->getObjValue());
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

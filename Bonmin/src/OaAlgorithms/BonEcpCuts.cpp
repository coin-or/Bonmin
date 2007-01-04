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
  double violation = nlp_->getConstraintViolation(
                     si.getColSolution(), si.getObjValue());
  //  std::cout<<"Constraint violation: "<<violation<<std::endl;
  if(violation <= 1e-02)
    return;
  solverManip * lpManip = NULL;
  bool infeasible = false;
  for(int i = 0 ; i < numRounds_ ; i++)
  {
    if( violation > 1e-02)
    {
      int numberCuts =  - cs.sizeRowCuts();
      nlp_->getOuterApproximation(cs, si.getColSolution(), 1);
      numberCuts += cs.sizeRowCuts();
      if(numberCuts > 0 && i + 1 < numRounds_){
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
        violation =  nlp_->getConstraintViolation(
                     lpManip->si()->getColSolution(), 
                     lpManip->si()->getObjValue());
	//        std::cout<<"Constraint violation: "<<violation<<std::endl;
         }
      else break;
    }
    else break;
  }
  if(!infeasible){
    lpManip->si()->resolve();
    if(lpManip->si()->isProvenPrimalInfeasible())
      objValue_ = 2e50;
    else
      objValue_ = lpManip->si()->getObjValue();}
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

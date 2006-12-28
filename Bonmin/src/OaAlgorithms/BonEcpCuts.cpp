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

void
EcpCuts::generateCuts(const OsiSolverInterface &si, OsiCuts & cs,
                      const CglTreeInfo info) const
{
  std::cout<<"Start ecp cut generation"<<std::endl;
  double violation = nlp_->getConstraintViolation(
                     si.getColSolution(), si.getObjValue());
  std::cout<<"Constraint violation: "<<violation<<std::endl;
  if(violation <= 1e-01)
    return;
  int numIt = 1; 
  solverManip * lpManip = NULL;
  for(int i = 0 ; i < numIt ; i++)
  {
    if( violation > 1e-01)
    {
      int numberCuts =  - cs.sizeRowCuts();
      nlp_->getOuterApproximation(cs, si.getColSolution(), 1);
      numberCuts += cs.sizeRowCuts();
      if(numberCuts > 0 && i + 1 < numIt){
        if(lpManip==NULL) lpManip = new solverManip(si);
        lpManip->installCuts(cs,numberCuts);
        lpManip->si()->resolve();
        violation =  nlp_->getConstraintViolation(
                     lpManip->si()->getColSolution(), 
                     lpManip->si()->getObjValue());
        std::cout<<"Constraint violation: "<<violation<<std::endl;
         }
      else
        break;
    }
    break;
  }
  delete lpManip;
  std::cout<<"End ecp cut generation"<<std::endl;
  return;
}

} // end namespace bonmin.

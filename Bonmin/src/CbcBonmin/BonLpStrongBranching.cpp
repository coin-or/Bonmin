#include "BonLpStrongBranching.hpp"
#include "BonEcpCuts.hpp"
#include "OsiClpSolverInterface.hpp"


namespace Bonmin{

  LpStrongBranching::LpStrongBranching(
          OsiTMINLPInterface * solver):
  BonChooseVariable(solver),
  maxCuttingPlaneIteration_(5){
  }
  
  LpStrongBranching::LpStrongBranching(
          const LpStrongBranching &other):
    BonChooseVariable(other),
    maxCuttingPlaneIteration_(other.maxCuttingPlaneIteration_){
   }

   LpStrongBranching&
   LpStrongBranching::operator=
                       (const LpStrongBranching& rhs){
     if(this != &rhs){
     BonChooseVariable::operator=(rhs);
     maxCuttingPlaneIteration_ = rhs. maxCuttingPlaneIteration_;
   }
   return *this;
   }

   OsiChooseVariable *
   LpStrongBranching::clone() const{
     return new LpStrongBranching(*this);
   }

  LpStrongBranching::~LpStrongBranching(){
  }
   int
   LpStrongBranching::fill_changes(OsiSolverInterface * solver,
                                   OsiBranchingInformation *info,
			       bool fixVariables, int numStrong,
			       double* change_down,
			       double* change_up, int& best_way)
   {
     std::cout<<"Start strong branching"<<std::endl;
     // get info about nlp solution and others
     const double * colsol = solver->getColSolution();
     const double * colLow = solver->getColLower();
     const double * colUp = solver->getColUpper();

     //Setup LP approximation model
     OsiTMINLPInterface * tminlpSi = dynamic_cast<OsiTMINLPInterface *>(solver);

     OsiClpSolverInterface lin;
     tminlpSi->extractLinearRelaxation(lin, 1, 0);
     lin.setDblParam(OsiDualObjectiveLimit, info->cutoff_);
     lin.resolve();
     CoinWarmStart * warm = lin.getWarmStart();
     double curObj = lin.getObjValue();
     EcpCuts ecp(tminlpSi, 5);
     int return_value = -1;
     for(int i = 0 ; i < numStrong ; i++){
       int & index = list_[i];
       const OsiObject * object = solver->object(index);
       int colnum = object->columnNumber();
       DBG_ASSERT(colnum != -1);

       double saveBound = colLow[colnum];
       double newBound = Min(ceil(colsol[colnum]), colUp[colnum]);
       lin.setColLower(colnum, newBound);
       lin.setWarmStart(warm);
       change_up[i] = ecp.doEcpRounds(lin, true) - curObj;
       if(change_up[i] >= 1e50){//Problem is infeasible force branch
         best_way = 1;
         return_value = i;
         break;
       }
       lin.setColLower(colnum, saveBound);
       
       newBound = Max(floor(colsol[colnum]), colLow[colnum]);
       lin.setColUpper(colnum, newBound);
       lin.setWarmStart(warm);
       change_down[i] = ecp.doEcpRounds(lin, true) - curObj;
       if(change_down[i] >= 1e50){//Problem is infeasible force branch
         best_way = 0;
         return_value = i;
         break;
       }
       lin.setColUpper(colnum, saveBound);
     } 
     std::cout<<"End strong branching"<<std::endl;
     delete warm;
     return return_value;
   }
}/* End namespace Bonmin. */

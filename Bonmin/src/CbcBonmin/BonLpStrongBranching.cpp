#include "BonLpStrongBranching.hpp"
#include "BonEcpCuts.hpp"
#include "OsiClpSolverInterface.hpp"


namespace Bonmin{

  LpStrongBranching::LpStrongBranching(
          OsiTMINLPInterface * solver):
  BonChooseVariable(solver),
  maxCuttingPlaneIterations_(5){
  }
  
  LpStrongBranching::LpStrongBranching(
          const LpStrongBranching &other):
    BonChooseVariable(other),
    maxCuttingPlaneIterations_(other.maxCuttingPlaneIterations_){
   }

   LpStrongBranching&
   LpStrongBranching::operator=
                       (const LpStrongBranching& rhs){
     if(this != &rhs){
     BonChooseVariable::operator=(rhs);
     maxCuttingPlaneIterations_ = rhs. maxCuttingPlaneIterations_;
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
     // get info about nlp solution and others
     const double * colsol = solver->getColSolution();
     const double * colLow = solver->getColLower();
     const double * colUp = solver->getColUpper();

     //Setup LP approximation model
     OsiTMINLPInterface * tminlpSi = dynamic_cast<OsiTMINLPInterface *>(solver);

     OsiClpSolverInterface lin;
     tminlpSi->extractLinearRelaxation(lin, 1, 0);
     lin.setDblParam(OsiDualObjectiveLimit, info->cutoff_);
     lin.messageHandler()->setLogLevel(0);
     lin.resolve();
     CoinWarmStart * warm = lin.getWarmStart();
     double curObj = lin.getObjValue();
     EcpCuts ecp(tminlpSi, maxCuttingPlaneIterations_);
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
       lin.resolve();
       if(lin.isProvenPrimalInfeasible())
	 {
	   best_way = 1;
	   return_value = i;
	   break;
	 }
       if(maxCuttingPlaneIterations_ > 0)
	 change_up[i] = ecp.doEcpRounds(lin, true) - curObj;
       else
	 change_up[i] = lin.getObjValue() - curObj;
       if(change_up[i] >= 1e50){//Problem is infeasible force branch
         best_way = 1;
         return_value = i;
         break;
       }
       lin.setColLower(colnum, saveBound);
       
       saveBound = colUp[colnum];
       newBound = Max(floor(colsol[colnum]), colLow[colnum]);
       lin.setColUpper(colnum, newBound);
       lin.setWarmStart(warm);
       lin.resolve();
       if(lin.isProvenPrimalInfeasible())
	 {
	   best_way = 0;
	   return_value = i;
	   break;
	 }
       else if(maxCuttingPlaneIterations_ > 0)
	 change_down[i] = ecp.doEcpRounds(lin, true) - curObj;
       else
	 change_down[i] = lin.getObjValue() - curObj;
       if(change_down[i] >= 1e50){//Problem is infeasible force branch
         best_way = 0;
         return_value = i;
         break;
       }
       lin.setColUpper(colnum, saveBound);
     }
     delete warm;
     return return_value;
   }
}/* End namespace Bonmin. */

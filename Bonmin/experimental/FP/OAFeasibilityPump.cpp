// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006 
// All Rights Reserved.
// This code is published under the Common Public License.
//
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 06/18/2005
// The ideas developped in this code are co-authored with: G. Cornuejols, A. Lodi, F. Margot and can be found 
// in:
// "A Feasibility Pump for Mixed Integer Nonlinear Programming" 
//  IBM Research Report RC 23682 February 2006
//

#include "BonminConfig.h"

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <iomanip>
#include <sstream>

//#define EXTRACT_LIN_REL

// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcBranchUser.hpp"
#include "CbcCompareUser.hpp"
#include "CbcCompareActual.hpp"
#include "CbcCutGenerator.hpp"
//#include "CbcHeuristicUser.hpp"
#include "BonAmplInterface.hpp"
#include "BonDummyHeuristic.hpp"
#include "BonOACutGenerator2.hpp"

//Use heuristics
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristicGreedy.hpp"

#include "BonAmplTMINLP.hpp"

#include "CglGomory.hpp"
//#include "CglProbing.hpp"

//#include "ClpQuadInterface.hpp"

// Heuristics would need adapting

//#include "CbcHeuristic.hpp"
//#define COIN_HAS_CPX
#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#else
#include "OsiCbcSolverInterface.hpp"
#endif

// Time
#include "CoinTime.hpp"
#include "FP.hpp"

using namespace Bonmin;

OptParam params;
static double BeginTimeGLOB;

/// solve problem with CPLEX for :
/// at least minNodes
/// then check every nodeInterval if a solution has been found and stop
/// do that for maxTime
#ifdef COIN_HAS_CPX
int findGoodSolution(OsiCpxSolverInterface &mip, double maxTime,
                     int& nTotalNodes)
{
  double beginTime = CoinCpuTime();
  int saveIntSolLim;
  int saveNodLim;
  double saveTiLim;
  CPXCENVptr env2 = mip.getEnvironmentPtr();
  CPXCLPptr lp = mip.getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL);
  CPXgetintparam(env2,CPX_PARAM_INTSOLLIM, &saveIntSolLim);
  CPXgetintparam(env2,CPX_PARAM_NODELIM, &saveNodLim);
  CPXgetdblparam(env2,CPX_PARAM_TILIM, &saveTiLim);

  CPXsetintparam(mip.getEnvironmentPtr(),CPX_PARAM_INTSOLLIM, 10000);
  CPXsetintparam(mip.getEnvironmentPtr(),CPX_PARAM_NODELIM, params.minNodes_);
  CPXsetdblparam(mip.getEnvironmentPtr(),CPX_PARAM_TILIM, maxTime);

  mip.branchAndBound();
  nTotalNodes += CPXgetnodecnt(env2,lp);
  int mipstat = CPXgetstat(env2,lp);
  CPXsetintparam(mip.getEnvironmentPtr(),CPX_PARAM_INTSOLLIM, 1);
  CPXsetintparam(mip.getEnvironmentPtr(),CPX_PARAM_NODELIM, params.nodeInterval_);
  while (mipstat == CPXMIP_NODE_LIM_INFEAS) {
    CPXsetdblparam(mip.getEnvironmentPtr(),CPX_PARAM_TILIM, maxTime - CoinCpuTime() + beginTime);
    mip.branchAndBound();
    nTotalNodes += CPXgetnodecnt(env2,lp);
    mipstat = CPXgetstat(env2,lp);
  }
  CPXsetintparam(mip.getEnvironmentPtr(),CPX_PARAM_INTSOLLIM, saveIntSolLim);
  CPXsetintparam(mip.getEnvironmentPtr(),CPX_PARAM_NODELIM, saveNodLim);
  CPXsetdblparam(mip.getEnvironmentPtr(),CPX_PARAM_TILIM, saveTiLim);
  return mipstat;
}
#else
int findGoodSolution(OsiSolverInterface &mip, int&nodeNumber)
{
  std::cerr<<"Does not work without CPX"<<std::endl;
  throw -1;
}

#endif
struct ResolutionInformation
{
  double time;
  double nlp_time;
  double mip_time;
  int n_iterations;
  int n_nlp_iterations;
  int n_mip_nodes;
  ResolutionInformation():time(0.), nlp_time(0.), mip_time(0.), n_iterations(0), n_nlp_iterations(0), n_mip_nodes(0)
  {}
  ResolutionInformation& operator +=(ResolutionInformation &other)
  {
    time +=other.time;
    nlp_time+=(other.nlp_time);
    mip_time+=(other.mip_time);
    n_iterations+=(other.n_iterations);
    n_nlp_iterations+=other.n_nlp_iterations;
    n_mip_nodes+=other.n_mip_nodes;
    return *this;
  }
};

#if 0
void
writeBoundFiles(AmplInterface& nlp, const double * originalLower, const double * originalUpper)
{
  const double * currentLower = nlp.getColLower();
  const double * currentUpper = nlp.getColUpper();

  CoinRelFltEq eq;
  std::string fBoundsName;
  nlp.getStrParam(OsiProbName,fBoundsName);
  fBoundsName+="_bounds_FP";
  std::string fModName = fBoundsName + ".mod";
  std::ofstream fBounds;
  std::ofstream fMod;
  bool hasVarNames = 0;
  if(nlp.getVarNames()!=NULL )
    hasVarNames=1;
  if(hasVarNames)
    fMod.open(fModName.c_str());
  fBounds.open(fBoundsName.c_str());

  for(int i = 0 ; i < nlp.getNumCols() ; i++) {

    if(nlp.isInteger(i) && !eq(currentLower[i],originalLower[i])) {
      if(hasVarNames)
        fMod<<"bounds"<<i<<": "
        <<nlp.getVarNames()[i]<<" >= "
        <<currentLower[i]<<";\n";
      fBounds<<"LO"<<"\t"<<i<<"\t"<<currentLower[i]<<std::endl;
    }
    if(nlp.isInteger(i) && !eq(currentUpper[i],originalUpper[i])) {
      if(hasVarNames)
        fMod<<"bounds"<<i<<": "
        <<nlp.getVarNames()[i]<<" <= "
        <<currentUpper[i]<<";\n";
      fBounds<<"UP"<<"\t"<<i<<"\t"<<currentUpper[i]<<std::endl;
    }
  }
}
#endif

double FP(AmplInterface &nlp,
          OsiSolverInterface &linearModel,
          int numIntCols, int * inds, double * vals,
          double maxTime, int maxIter, ResolutionInformation& info,
          double ub, int &provenInfeas, double *&solution)
{
  provenInfeas=0;
#ifdef COIN_HAS_CPX

  OsiCpxSolverInterface * cpxSi = dynamic_cast<OsiCpxSolverInterface *>
                                  (&linearModel);
  OsiCpxSolverInterface &mip = *cpxSi;
#endif

  int &nMajorIt = info.n_iterations;
  int &nTotalNodes = info.n_mip_nodes;
  int &nNlpIterations = info.n_nlp_iterations;
  double& bbTime = info.mip_time;
  double& nlpTime = info.nlp_time;
  info.time = -CoinCpuTime();
  double objValue = DBL_MAX;
  {
    OsiCuts cs;
    //First get the closest point to the current integer (NLP-infeasible) optimum
    nlpTime -= CoinCpuTime();
    double dist = nlp.getFeasibilityOuterApproximation( numIntCols, vals, inds, cs, false, true);
    nlpTime += CoinCpuTime();
    if(dist < 1e-06)//do integer infeasibility check on variables
    {
      //Something wrong shouldn't be
      std::cout<<"Feasibility subproblem has objective 0 while problem was claimed infeasible before"<<std::endl
      <<"I am confused, exiting with error"<<std::endl;
      throw -1;
    }
    OsiSolverInterface::ApplyCutsReturnCode acRc;

    //       linearModel.writeMps("test1");
    //       acRc = linearModel.applyCuts(cs);
    //       linearModel.writeMps("test2");

    // Print applyCuts return code
    std::cout <<cs.sizeCuts() <<" cuts were generated" <<std::endl;
    std::cout <<"  " <<acRc.getNumInconsistent() <<" were inconsistent" <<std::endl;
    std::cout <<"  " <<acRc.getNumInconsistentWrtIntegerModel()
    <<" were inconsistent for this problem" <<std::endl;
    std::cout <<"  " <<acRc.getNumInfeasible() <<" were infeasible" <<std::endl;
    std::cout <<"  " <<acRc.getNumIneffective() <<" were ineffective" <<std::endl;
    std::cout <<"  " <<acRc.getNumApplied() <<" were applied" <<std::endl;
    std::cout << std::endl << std::endl;
  }
  //Modify the linearModel for feasibility pump
  //save objective
  int numcols = linearModel.getNumCols();
  double * saveObj = new double[numcols];
  CoinCopyN(linearModel.getObjCoefficients(), numcols, saveObj);
  //double * obj = new double[2*numIntCols];
  //int k = 0;

  //if ub is given add a "cutoff" row
  if(ub < 1e100) {
    //change bound on alpha
    linearModel.setColUpper(numcols-1,ub);
    //     addCutOff=1;
    //     CoinPackedVector * rowToAdd = (CoinPackedVector *) emptyCols[k];
    //     rowToAdd->insert(origNumCols-1,1.);
    //     colUb[k] = ub;
    //     colLb[k] = -linearModel.getInfinity();
  }
  for(int i = 0;i< numcols ;i++)
    if(!linearModel.isInteger(i))
      linearModel.setObjCoeff(i,0.);

  //done
  int nAddedCuts = 0;

  bool solved = 0;
  bool numericFailure = 0;
  int numScaleIncrease = 0;
  while((nMajorIt < maxIter) && (!solved && ( CoinCpuTime() + info.time < maxTime) && !numericFailure && numScaleIncrease < 5)) {
    //Change the objective of the MIP to get the closest point to rounding of NLP optimum
    for(int i = 0; i < numIntCols; i++) {
      linearModel.setObjCoeff(inds[i],1 - 2* nlp.getColSolution()[inds[i]]);
    }
#ifdef COIN_HAS_CPX
    //	      CPXsetintparam(cpxSi->getEnvironmentPtr(),CPX_PARAM_NODELIM, nMaxNodes - nTotalNodes);
    //      CPXsetintparam(cpxSi->getEnvironmentPtr(),CPX_PARAM_INTSOLLIM, 1);
    //      CPXsetintparam(mip.getEnvironmentPtr(),CPX_PARAM_HEURFREQ, 10);
    //      CPXsetdblparam(cpxSi->getEnvironmentPtr(),CPX_PARAM_TILIM, maxTime - CoinCpuTime() + time);

#else
    linearModel.initialSolve();
    CbcModel mip(linearModel);
    //        mip.writeMps("oa");
    CbcStrategyDefault defaultStrategy;
    mip.setStrategy(defaultStrategy);
    mip.solver()->messageHandler()->setLogLevel(0);

    //Add some heuristics to get feasible solutions
    mip.setMaximumSeconds(maxTime - CoinCpuTime() - info.time);
    mip.setMaximumSolutions(1);
#endif

    bbTime -= CoinCpuTime();
#ifdef COIN_HAS_CPX

    int mipstat = findGoodSolution(mip, maxTime - CoinCpuTime() - info.time, nTotalNodes);
#else
    //  mip.branchAndBound();
#endif

    bbTime += CoinCpuTime();
    //      mip.writeMps("oa1");
    OsiCuts cs;
    nMajorIt++;

#ifdef COIN_HAS_CPX

    const double * colsol = mip.getColSolution();
    //int nNodes = 0;
    nTotalNodes += CPXgetnodecnt(cpxSi->getEnvironmentPtr(),cpxSi->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
    //    int mipstat = CPXgetstat(mip.getEnvironmentPtr(),mip.getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
    //   CPXsetintparam(mip.getEnvironmentPtr(),CPX_PARAM_HEURFREQ, 0);
    if(mipstat == CPXMIP_INFEASIBLE) {
      provenInfeas=1;
      info.time += CoinCpuTime();
      for(int i = 0 ; i < numcols; i++)
        linearModel.setObjCoeff(i, saveObj[i]);
      return 1e75;
    }
    else if(mipstat == CPXMIP_TIME_LIM_FEAS || mipstat == CPXMIP_TIME_LIM_INFEAS) {
      provenInfeas= -1;
      info.time += CoinCpuTime();
      for(int i = 0 ; i < numcols; i++)
        linearModel.setObjCoeff(i, saveObj[i]);
      return 1e75;
    }
#else
    if(mip.bestSolution() == NULL)
      break;
    if(mip.isNodeLimitReached()) {
      break;
    }
    //      if(!mip.isSolutionLimitReached())
    //          break;
    nTotalNodes += mip.getNodeCount();
    const double * colsol = mip.bestSolution();
#endif

    for(int i = 0 ; i < numIntCols ; i++) {
      vals[i] = floor(colsol[inds[i]] + 0.5);
      vals[i] = max(mip.getColLower()[inds[i]],vals[i]);
      vals[i] = min(mip.getColUpper()[inds[i]],vals[i]);
    }
    nlpTime -= CoinCpuTime();
    double dist = nlp.getFeasibilityOuterApproximation( numIntCols, vals, inds, cs, false, true);
    nlpTime += CoinCpuTime();
    nNlpIterations += nlp.getIterationCount();
    if(!nlp.isProvenOptimal()) {
      std::cout<<"?????"<<std::endl;
      throw -1;
    }
    nlpTime += CoinCpuTime();
    if(dist < 1e-04)//do integer infeasibility check on variables
    {
      solved = 1;
      double norm_inf = DBL_MAX;
      for(int i = 0 ; i < nlp.getNumCols() && solved; i++)
      {
        if(nlp.isInteger(i)) {
          const double &value = nlp.getColSolution()[i];
          double IIf = fabs( value- floor(value + 0.5));
          norm_inf = min(norm_inf, IIf);
          if(fabs( value- floor(value + 0.5)) > 1e-06) {
            std::cout<<"Variable "<<i<<" has integer infeasibility : "<<IIf<<std::endl;
            solved=0;
            numericFailure = 1;
#if 0

            numScaleIncrease++;
            nlp.fpnlp()->setObjectiveScaling(10* nlp.fpnlp()->getObjectiveScaling());
            for(int i = 0 ; i < numIntCols ; i++) {
              nlp.setColLower(inds[i], linearModel.getColLower()[inds[i]]);
              nlp.setColUpper(inds[i], linearModel.getColUpper()[inds[i]]);
            }
#endif

          }
        }
      }
      std::cout<<"Found a solution with maximal integer infeasibility of "<<norm_inf<<std::endl;
    }

    OsiSolverInterface::ApplyCutsReturnCode acRc;
    std::cout<<"Cut violation :"<<cs.rowCut(0).violated(colsol)
    <<std::endl;
    //       cs.printCuts();

    linearModel.writeMps("test1");
    acRc = linearModel.applyCuts(cs);
    linearModel.writeMps("test2");
    nAddedCuts += cs.sizeRowCuts();
    // Print applyCuts return code
    std::cout <<cs.sizeCuts() <<" cuts were generated" <<std::endl;
    std::cout <<"  " <<acRc.getNumInconsistent() <<" were inconsistent" <<std::endl;
    std::cout <<"  " <<acRc.getNumInconsistentWrtIntegerModel()
    <<" were inconsistent for this problem" <<std::endl;
    std::cout <<"  " <<acRc.getNumInfeasible() <<" were infeasible" <<std::endl;
    std::cout <<"  " <<acRc.getNumIneffective() <<" were ineffective" <<std::endl;
    std::cout <<"  " <<acRc.getNumApplied() <<" were applied" <<std::endl;
    std::cout << std::endl << std::endl;


    if(solved) {
      //Set warm start point to the last point found (which is feasible for this relaxation)
      nlp.setColSolution(nlp.getColSolution());
      nlp.setRowPrice(nlp.getRowPrice());
      nlp.solver()->enableWarmStart();
      //Resolve the NLP with fixed variables and original objective function
      for(int i = 0; i < numIntCols; i++) {
        nlp.setColLower(inds[i], vals[i]);
        nlp.setColUpper(inds[i], vals[i]);
      }
      //nlp.turnOnSolverOutput();
      nlp.initialSolve();
      if(nlp.isProvenOptimal()) {
        OsiCuts cs;
        nlp.getOuterApproximation(cs,1, NULL, true);
        linearModel.applyCuts(cs);
        nAddedCuts += cs.sizeRowCuts();
        std::cout<<" FP found easible solution of value "<<nlp.getObjValue()<<" in "<< CoinCpuTime() - BeginTimeGLOB<<" seconds, "<<nMajorIt<<" major iterations, took"
        <<nTotalNodes<<" nodes."<<std::endl;
        objValue = nlp.getObjValue();
        if (solution==NULL) solution = new double[numcols];
        CoinCopyN(nlp.getColSolution(), numcols, solution);
      }
      else {
        std::cout<<" FP converged in "<<CoinCpuTime()- BeginTimeGLOB<<" seconds, "<<nMajorIt<<" major iterations, took"
        <<nTotalNodes<<" nodes, but solution seems infeasible"<<std::endl;
        std::cout<<"Increasing scaling of objective and restarting"<<std::endl;
        solved = 0;
        numScaleIncrease++;
        numericFailure = 1;
#if 0

        nlp.fpnlp()->setObjectiveScaling(10* nlp.fpnlp()->getObjectiveScaling());
        for(int i = 0 ; i < numIntCols ; i++) {
          nlp.setColLower(inds[i], linearModel.getColLower()[inds[i]]);
          nlp.setColUpper(inds[i], linearModel.getColUpper()[inds[i]]);
        }
#endif

      }

    }

  }

  if(!solved) {
    if(numericFailure)
      std::cout<<"FP aborted because of a numberical failure : "<<CoinCpuTime() - BeginTimeGLOB<<",  "<<nMajorIt<<" major iterations, took"
      <<nTotalNodes<<" nodes."<<std::endl;
    else
      std::cout<<"FP aborted on time limit elapsed time : "<<CoinCpuTime() - BeginTimeGLOB<<",  "<<nMajorIt<<" major iterations, took"
      <<nTotalNodes<<" nodes."<<std::endl;
    provenInfeas = -1;

  }
  std::cout<<"Nlp solve time : "<<nlpTime<<" B&B solve time : "<<bbTime<<std::endl;
  //    restore objective
  for(int i = 0 ; i < numcols; i++)
    linearModel.setObjCoeff(i, saveObj[i]);
  delete [] saveObj;
  info.time += CoinCpuTime();
  return objValue;
}


#if 0 // Old enhanced OA code without the iFP at beginning
int main3 (int argc, char *argv[])
{
  std::cout.precision(11);
  using namespace Ipopt;

  // Read in model using argv[1]
  char * pbName = new char[strlen(argv[1])+1];
  strcpy(pbName, argv[1]);
  BonminAmplInterface solver1(argv);
  solver1.turnOnSolverOutput();
  bool doFp=true;

#ifdef COIN_HAS_CPX

  OsiCpxSolverInterface solver;
  OsiCpxSolverInterface &mip = solver;
#else

  OsiClpSolverInterface solver;
#endif

  solver.messageHandler()->setLogLevel(0);


  //Setup timers since an nlp is solved to extract the linear relaxation
  BeginTimeGLOB= CoinCpuTime();
  double bbTime = 0;
  double nlpTime = - BeginTimeGLOB;

  ResolutionInformation FP_infos;
  solver1.extractLinearRelaxation(solver, 1);
  nlpTime += CoinCpuTime();

#ifdef EXTRACT_LIN_REL

  std::string linearName="Lin";
  linearName += pbName;
  solver.writeMpsNative(linearName.c_str(),NULL,NULL);
#endif

  int numIntCols = 0;
  int * inds = new int[solver1.getNumCols()]; //indices of integer variables
  double * x = new double[solver1.getNumCols()]; //to store values of integer variables later on
  //int origNumCols = solver.getNumCols();
  //int origNumRows = solver.getNumRows();
  for(int i = 0; i < solver1.getNumCols(); i++) {
    if(solver1.isInteger(i)) {
      inds[numIntCols++] = i;
    }
  }

  double lb = solver1.getObjValue();
  double ub = DBL_MAX;
  //Limits of the procedure
  //  int nMaxNodes = 100000;
  params.maxTime_ = 2.*3600.;

  //counters for iterations, nodes,...
  int nMajorIt = 0;
  int nTotalNodes = 0;
  int nTimesFPCalled =0;
  double precision = 1e-02;
  //bool solved = 0;
  double firstIterationTime;
  solver.messageHandler()->setLogLevel(0);
  int numNotFound = 0;
  bool minutePassed = 0;
  while( (ub - lb) > precision && (CoinCpuTime() - BeginTimeGLOB < params.maxTime_) ) {
    solver.resolve();
#ifndef COIN_HAS_CPX

    CbcModel mip(solver);
    //  solver.writeMps("oa");
    CbcStrategyDefault defaultStrategy(1,8,4);
    mip.setStrategy(defaultStrategy);

    mip.setMaximumSeconds(params.maxTime_ - CoinCpuTime() + BeginTimeGLOB);
    mip.setMaximumSolutions(3);
    mip.setCutoff(ub);
    mip.setLogLevel(1);
    mip.solver()->messageHandler()->setLogLevel(0);
    //      mip.solver()->writeMps("test");
#else
    //	      CPXsetintparam(cpxSi->getEnvironmentPtr(),CPX_PARAM_NODELIM, nMaxNodes - nTotalNodes);
    //   CPXsetintparam(mip.getEnvironmentPtr(),CPX_PARAM_INTSOLLIM, 1);
    CPXsetdblparam(mip.getEnvironmentPtr(),CPX_PARAM_TILIM, params.maxTime_ - CoinCpuTime() + BeginTimeGLOB);
    CPXsetdblparam(mip.getEnvironmentPtr(),CPX_PARAM_CUTUP, ub);

#endif

    bbTime -= CoinCpuTime();

    //      mip.solver()->writeMps("oa");

    //int mipstat = findGoodSolution(mip, minNodes,
    //                             nodeInterval, maxTime - CoinCpuTime() + time, nTotalNodes)

    mip.branchAndBound();
    bbTime += CoinCpuTime();
    if(nMajorIt==0) {
      firstIterationTime = bbTime;
    }
#ifdef EXTRACT_LIN_REL
    std::string oa_name=linearName;
    std::ostringstream t;
    t<<nMajorIt;
    oa_name+=t.str();
#endif

    OsiCuts cs;
#ifdef COIN_HAS_CPX

    const double *colsol=mip.getColSolution();
    nTotalNodes += CPXgetnodecnt(mip.getEnvironmentPtr(),mip.getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
    double new_lb = 0;
    int mipstat = CPXgetstat(mip.getEnvironmentPtr(),mip.getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
    //    if(mipstat ==
    if((mipstat != CPXMIP_OPTIMAL && mipstat != CPXMIP_OPTIMAL_TOL)) {
      int status = CPXgetbestobjval(mip.getEnvironmentPtr(),mip.getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL), &new_lb);
      lb=max(lb,new_lb);
      if(status)
        throw CoinError("Error in getting CPLEX best bound","IpCbcOACutGenerator2","siBestObj");
    }
    else
      lb = mip.getObjValue();

    CPXsetdblparam(mip.getEnvironmentPtr(),CPX_PARAM_CUTUP, mip.getInfinity());

#else
    //      if(!mip.isSolutionLimitReached())
    //          break;
    nTotalNodes += mip.getNodeCount();
    const double * colsol = mip.bestSolution();
    lb = mip.getBestPossibleObjValue();
    if(mip.bestSolution() == NULL)
      break;
#endif

    nMajorIt++;

    std::cout<<"Found new lower bound of : "<<lb<<std::endl;
    if ((ub - lb) < precision)
      break;
    for(int i = 0 ; i < numIntCols ; i++) {
      x[i] = max(solver.getColLower()[inds[i]],colsol[inds[i]]);
      x[i] = min(solver.getColUpper()[inds[i]],x[i]);
      solver1.setColLower(inds[i], x[i]);
      solver1.setColUpper(inds[i], x[i]);
    }
    solver1.turnOnSolverOutput();
    solver1.initialSolve();
    if(solver1.isProvenOptimal()) {
      std::cout<<pbName<<" OA found easible solution of value "
      <<solver1.getObjValue()<<" in "
      <<CoinCpuTime() -BeginTimeGLOB<<" seconds, "
      <<nMajorIt<<" major iterations, took"
      <<nTotalNodes<<" nodes."<<std::endl;
      if(ub > 1e100)
        std::cout<<"First solution "<< solver1.getObjValue()<<" after "
        <<CoinCpuTime() -BeginTimeGLOB<<" secs."<<std::endl;
      if(!minutePassed && (CoinCpuTime() -BeginTimeGLOB) > 61.) {
        minutePassed = true;
        std::cout<<"Solution at minute "<< ub<<" after "
        <<CoinCpuTime() -BeginTimeGLOB<<" secs."<<std::endl;
      }

      OsiCuts cs;
      solver1.getOuterApproximation(cs,1);
      OsiSolverInterface::ApplyCutsReturnCode acRc;
      ub = min (solver1.getObjValue(), ub);
      acRc = solver.applyCuts(cs);
      // Print applyCuts return code
      std::cout <<cs.sizeCuts() <<" cuts were generated" <<std::endl;
      std::cout <<"  " <<acRc.getNumInconsistent() <<" were inconsistent" <<std::endl;
      std::cout <<"  " <<acRc.getNumInconsistentWrtIntegerModel()
      <<" were inconsistent for this problem" <<std::endl;
      std::cout <<"  " <<acRc.getNumInfeasible() <<" were infeasible" <<std::endl;
      std::cout <<"  " <<acRc.getNumIneffective() <<" were ineffective" <<std::endl;
      std::cout <<"  " <<acRc.getNumApplied() <<" were applied" <<std::endl;
      std::cout << std::endl << std::endl;

    }
    else if(solver1.isAbandoned() || solver1.isIterationLimitReached()) {
      TNLPSolver::UnsolvedError * E=solver1.newUnsolvedError(0,argv[1],"not solved");
      E->writeDiffFiles();
      solver1.turnOnSolverOutput();
      solver1.initialSolve();

      std::cerr<<"Error"<<std::endl;
      throw -1;
    }

    else {
      if (doFp)//Launch FP
      {
        //Reset NLP bounds
        for(int i = 0 ; i < numIntCols ; i++)
        {
          solver1.setColLower(inds[i], solver.getColLower()[inds[i]]);
          solver1.setColUpper(inds[i], solver.getColUpper()[inds[i]]);
        }
        ResolutionInformation infos;
        nTimesFPCalled++;
        int provenInfeas =0;

        //compute FP maxTime
        //	      firstIterationTime = 1.;
        double fpTime = min( params.maxTime_ - CoinCpuTime() + BeginTimeGLOB,120.);
        int fpMaxIter = 30;
        ub = min(ub, FP(solver1, solver,numIntCols, inds, x,fpTime, fpMaxIter,  infos, ub*(1-1e-03), provenInfeas, solution));
        if(provenInfeas == 1)
          lb = mip.getInfinity();
        if(provenInfeas == -1)
        {
          numNotFound++;
          if(numNotFound==2)
            doFp = 0;
        }
        // 	      if(infos.n_iterations >= 10)
        // 		doFp = 0;
        FP_infos+=infos;
      }
      else {
        std::cout<<"Adding feasibility cuts based on 1-norm of"
        <<"constraint satisfaction"<<std::endl;
        // solver1.getOuterApproximation(cs,1);
        solver1.getFeasibilityOuterApproximation( numIntCols, x, inds, cs);
        OsiSolverInterface::ApplyCutsReturnCode acRc;
        acRc = solver.applyCuts(cs);
        // Print applyCuts return code
        std::cout <<cs.sizeCuts() <<" cuts were generated" <<std::endl;
        std::cout <<"  " <<acRc.getNumInconsistent() <<" were inconsistent" <<std::endl;
        std::cout <<"  " <<acRc.getNumInconsistentWrtIntegerModel()
        <<" were inconsistent for this problem" <<std::endl;
        std::cout <<"  " <<acRc.getNumInfeasible() <<" were infeasible" <<std::endl;
        std::cout <<"  " <<acRc.getNumIneffective() <<" were ineffective" <<std::endl;
        std::cout <<"  " <<acRc.getNumApplied() <<" were applied" <<std::endl;
        std::cout << std::endl << std::endl;
      }
    }

  }
  double time = CoinCpuTime() - BeginTimeGLOB;
  std::cout<<pbName<<" FP driven OA found optimal solution of value "<<ub<<" (lower bound is "<<lb<<") in "<<time<<" seconds, "<<nMajorIt<<" major iterations, "
  <<nTimesFPCalled<<" calls to FP ( "
  <<FP_infos.n_iterations <<" it, and "<< FP_infos.time
  <<"sec) overall took "
  <<nTotalNodes + FP_infos.n_mip_nodes<<" nodes."<<std::endl;
  std::cout<<"Nlp solve time : "<<nlpTime + FP_infos.nlp_time<<" B&B solve time : "<<bbTime + FP_infos.mip_time<<std::endl;
  delete[] inds;
  delete[] x;
  delete [] pbName;
  return 0;
}
#endif




/** Enhanced OA code */
int enhancedOA(AmplInterface & solver1, bool doFp,
               double *& solution)
{
  bool nonConvex = 0;
  std::string pbName;
  solver1.getStrParam(OsiProbName, pbName);
#ifdef COIN_HAS_CPX

  OsiCpxSolverInterface solver;
  OsiCpxSolverInterface &mip = solver;
#else

  OsiClpSolverInterface solver;
#endif

  solver.messageHandler()->setLogLevel(0);


  //Setup timers since an nlp is solved to extract the linear relaxation
  BeginTimeGLOB= CoinCpuTime();
  double bbTime = 0;
  double nlpTime = - BeginTimeGLOB;

  ResolutionInformation FP_infos;
  solver1.extractLinearRelaxation(solver, !nonConvex);
  //solver1.fpnlp()->setObjectiveScaling(1000* solver1.fpnlp()->getObjectiveScaling());
  nlpTime += CoinCpuTime();
  int numcols = solver.getNumCols();
  double * saveObj = new double[numcols];
  CoinCopyN(solver.getObjCoefficients(), numcols, saveObj);

#ifdef EXTRACT_LIN_REL
  {
    std::string linearName="Lin";
    linearName += basename(pbName);
    solver.writeMpsNative(linearName.c_str(),NULL,NULL);
  }
#endif

  int numIntCols = 0;
  int * inds = new int[solver1.getNumCols()]; //indices of integer variables
  double * x = new double[solver1.getNumCols()]; //to store values of integer variables later on
  //int origNumCols = solver.getNumCols();
  //int origNumRows = solver.getNumRows();
  for(int i = 0; i < solver1.getNumCols(); i++) {
    if(solver1.isInteger(i)) {
      inds[numIntCols++] = i;
      x[i]=0.5;
    }
  }

  double lb = solver1.getObjValue();
  double ub = DBL_MAX;
  //Limits of the procedure
  //  int nMaxNodes = 100000;
  double maxFpTime = 60.;//params.maxTime_+1;


  //counters for iterations, nodes,...
  int nMajorIt = 0;
  int nTotalNodes = 0;
  int nTimesFPCalled =0;
  int nTimesInitFPCalled =0;
  int nFPIterations =0;
  double precision = 1e-04;
  bool solved = 0;
  double firstIterationTime;
  solver.messageHandler()->setLogLevel(0);
  int numNotFound = 0;

  bool feasible = 1;

  double FPTime = 0.;
  
  if(doFp)
  {
  
  while(ub-lb> precision && feasible && ( CoinCpuTime() - BeginTimeGLOB < maxFpTime)) {
    nTimesInitFPCalled++;
    int nIt = 0;
    solved = 0;
    while(feasible & !solved &&
          (CoinCpuTime() - BeginTimeGLOB < maxFpTime)
          //   && nIt < 20
         ) {
      nIt++;
      nFPIterations++;
      for(int i = 0 ; i < numIntCols; i++) {
        solver.setObjCoeff(inds[i],1 - 2* solver1.getColSolution()[inds[i]]);
      }
      //change bound on alpha
      if(!nonConvex) {
        int sign = ub > 0 ? -1 : 1;
        solver.setColUpper(solver.getNumCols()-1,ub*(1 + sign *1e-03));
      }
#ifndef COIN_HAS_CPX
      solver.initialSolve();
      solver.messageHandler()->setLogLevel(0);

      CbcStrategyDefault defaultStrategy;
      CbcModel mip(solver);
      mip.setStrategy(defaultStrategy);
      mip.solver()->messageHandler()->setLogLevel(0);
      mip.setLogLevel(0);
      mip.setMaximumSeconds(params.maxTime_ - CoinCpuTime() - BeginTimeGLOB);
      mip.setMaximumSolutions(1);
#else
      //	      CPXsetintparam(cpxSi->getEnvironmentPtr(),CPX_PARAM_NODELIM, nMaxNodes - nTotalNodes);
    //  CPXsetintparam(mip.getEnvironmentPtr(),CPX_PARAM_INTSOLLIM, 3);
      CPXsetdblparam(mip.getEnvironmentPtr(),CPX_PARAM_TILIM,
                     maxFpTime - CoinCpuTime() + BeginTimeGLOB);
      
#endif

      bbTime -= CoinCpuTime();

#ifdef COIN_HAS_CPX

      int mipstat = findGoodSolution(mip, maxFpTime - CoinCpuTime() - BeginTimeGLOB, nTotalNodes);
#else

      mip.branchAndBound();
#endif

      bbTime += CoinCpuTime();
      //      mip.solver()->writeMps("oa");

      OsiCuts cs;
#ifdef COIN_HAS_CPX

      const double *colsol=mip.getColSolution();
      //	  nTotalNodes += CPXgetnodecnt(mip.getEnvironmentPtr(),mip.getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
      //int mipstat = CPXgetstat(mip.getEnvironmentPtr(),mip.getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
      if(mipstat == CPXMIP_INFEASIBLE || mipstat == CPXMIP_INForUNBD) {
        std::cout<<"I found an Infeasible mip"<<std::endl;
        feasible = 0;
        break;
      }

#else
      if(mip.bestSolution() == NULL)
        break;
      if(mip.isNodeLimitReached()) {
        break;
      }
      //		if(!mip.isSolutionLimitReached())
      //			break;
      nTotalNodes += mip.getNodeCount();
      const double * colsol = mip.bestSolution();
#endif

      nMajorIt++;
      for(int i = 0 ; i < numIntCols ; i++) {
        x[i] = max(solver1.getColLower()[inds[i]],colsol[inds[i]]);
        x[i] = min(solver1.getColUpper()[inds[i]],x[i]);
        //				std::cout<<"Var "<<ind[i]<<" value "<<x[i]<<std::endl;
      }
      nlpTime -= CoinCpuTime();
      double dist = solver1.getFeasibilityOuterApproximation( numIntCols, x, inds, cs, false, true);
      nlpTime += CoinCpuTime();
      std::cout<<"Dist : "<<dist<<std::endl;
      if(dist < 1e-04)//do integer infeasibility check on variables
      {
        solved = 1;
        for(int i = 0 ; i < solver1.getNumCols() && solved; i++)
        {
          if(solver1.isInteger(i)) {
            const double &value = solver1.getColSolution()[i];
            if(fabs( value- floor(value + 0.5)) > 1e-04) {
              solved = 0;
              std::cout<<"Variable "<<i<<" has integer infeasibility : "<<fabs( value- floor(value + 0.5))<<std::endl;
            }
          }
        }
        //	    feasible = !solved;
      }

      OsiSolverInterface::ApplyCutsReturnCode acRc;
      acRc = solver.applyCuts(cs);
      // Print applyCuts return code
      std::cout <<cs.sizeCuts() <<" cuts were generated" <<std::endl;
      std::cout <<"  " <<acRc.getNumInconsistent() <<" were inconsistent" <<std::endl;
      std::cout <<"  " <<acRc.getNumInconsistentWrtIntegerModel()
      <<" were inconsistent for this problem" <<std::endl;
      std::cout <<"  " <<acRc.getNumInfeasible() <<" were infeasible" <<std::endl;
      std::cout <<"  " <<acRc.getNumIneffective() <<" were ineffective" <<std::endl;
      std::cout <<"  " <<acRc.getNumApplied() <<" were applied" <<std::endl;
      std::cout << std::endl << std::endl;
    }
    double time = CoinCpuTime() - BeginTimeGLOB;
    if(solved) {
      //Resolve the NLP with fixed variables and original objective function
      for(int i = 0; i < numIntCols; i++) {
        solver1.setColLower(inds[i], x[i]);
        solver1.setColUpper(inds[i], x[i]);
      }
      solver1.turnOnSolverOutput();
      solver1.initialSolve();
      if(solver1.isProvenOptimal()) {
        ub = min(solver1.getObjValue(), ub);
        if (solution==NULL) solution = new double[numcols];
        CoinCopyN(solver1.getColSolution(), numcols, solution);
        std::cout<<pbName<<" FP found easible solution of value "
        <<solver1.getObjValue()<<" in "
        <<CoinCpuTime() - BeginTimeGLOB<<" seconds, "
        <<nMajorIt<<" major iterations, took"
        <<nTotalNodes<<" nodes."<<std::endl;
        if(nTimesInitFPCalled == 1)
          std::cout<<"FIRST sol found in : "<<time<<" seconds value is "<<solver1.getObjValue()
          <<" iterations "<<nFPIterations<<std::endl;

        solved = 1;
        OsiCuts cs;
        solver1.getOuterApproximation(cs,1, NULL, true);
        OsiSolverInterface::ApplyCutsReturnCode acRc;
        acRc = solver.applyCuts(cs);
        for(int i = 0 ; i < numIntCols ; i++) {
          solver1.setColLower(inds[i], solver.getColLower()[inds[i]]);
          solver1.setColUpper(inds[i], solver.getColUpper()[inds[i]]);
        }
        solver1.initialSolve();
      }
      else// may have cq problems add a feasibility cut
      {
        std::cout<<pbName<<" FP converged in "<<time<<" seconds, "<<nMajorIt<<" major iterations, took"
        <<nTotalNodes<<" nodes, but solution seems infeasible"<<std::endl;
        CoinPackedVector v;
        double rhs = 1.;
        for(int i = 0 ; i < numIntCols ; i++) {
          v.insert(inds[i], -(2*x[i] - 1));
          rhs -= x[i];
        }
        OsiCuts cs;
        OsiRowCut cut;
        cut.setRow(v);
        cut.setLb(rhs);
        cut.setUb(1e300);
        cut.setGloballyValid();
        cs.insert(cut);

        OsiSolverInterface::ApplyCutsReturnCode acRc;
        acRc = solver.applyCuts(cs);
        // Print applyCuts return code
        std::cout <<cs.sizeCuts() <<" cuts were generated" <<std::endl;
        std::cout <<"  " <<acRc.getNumInconsistent() <<" were inconsistent" <<std::endl;
        std::cout <<"  " <<acRc.getNumInconsistentWrtIntegerModel()
        <<" were inconsistent for this problem" <<std::endl;
        std::cout <<"  " <<acRc.getNumInfeasible() <<" were infeasible" <<std::endl;
        std::cout <<"  " <<acRc.getNumIneffective() <<" were ineffective" <<std::endl;
        std::cout <<"  " <<acRc.getNumApplied() <<" were applied" <<std::endl;
        std::cout << std::endl << std::endl;

      }
    }
    else if (feasible) {
      //restart from NLP-relaxation optimum
      solver1.initialSolve();
      nIt = 0;
    }
    if(nonConvex)
      break;
  }

  //revert mip to original objective function
  for(int i = 0 ; i < numcols; i++)
    solver.setObjCoeff(i, saveObj[i]);

  solver.setColUpper(solver.getNumCols()-1,ub);

  FPTime = CoinCpuTime() - BeginTimeGLOB;
  std::cout<<"Best known sol after : "<<FPTime<<" is "<<ub<<std::endl;

#ifdef EXTRACT_LIN_REL

  {
    std::string linearName="Lin1Min";
    linearName += basename(pbName);
    solver.writeMpsNative(linearName.c_str(),NULL,NULL);
  }
#endif

  solver.setColUpper(solver.getNumCols()-1,ub);
  
}
  while( (ub - lb) > precision && (CoinCpuTime() - BeginTimeGLOB < params.maxTime_) ) {
#ifdef EXTRACT_LIN_REL
    {
      if(CoinCpuTime() - BeginTimeGLOB > 5. * 60.)
      {
        std::string linearName="Lin5min";
        linearName += basename(pbName);
        solver.writeMpsNative(linearName.c_str(),NULL,NULL);
        return 0;
      }
    }
#endif

    solver.resolve();
#ifndef COIN_HAS_CPX

    solver.messageHandler()->setLogLevel(0);
    CbcModel mip(solver);
    //  solver.writeMps("oa");
    CbcStrategyDefault defaultStrategy(1,8,4);
    mip.setStrategy(defaultStrategy);

    mip.setMaximumSeconds(params.maxTime_ - CoinCpuTime() + BeginTimeGLOB);
    mip.setMaximumSolutions(3);
    mip.setCutoff(ub);
    mip.setLogLevel(1);
    mip.solver()->messageHandler()->setLogLevel(0);
    //      mip.solver()->writeMps("test");
#else
    //	      CPXsetintparam(cpxSi->getEnvironmentPtr(),CPX_PARAM_NODELIM, nMaxNodes - nTotalNodes);
    //   CPXsetintparam(mip.getEnvironmentPtr(),CPX_PARAM_INTSOLLIM, 1);
    CPXsetdblparam(mip.getEnvironmentPtr(),CPX_PARAM_TILIM, params.maxTime_ - CoinCpuTime() + BeginTimeGLOB);
    CPXsetdblparam(mip.getEnvironmentPtr(),CPX_PARAM_CUTUP, ub);

#endif

    bbTime -= CoinCpuTime();

    //      mip.solver()->writeMps("oa");
    //int mipstat = findGoodSolution(mip, minNodes,
    //                             nodeInterval, maxTime - CoinCpuTime() + time, nTotalNodes)

    mip.branchAndBound();
    bbTime += CoinCpuTime();
    if(nMajorIt==0) {
      firstIterationTime = bbTime;
    }

    OsiCuts cs;
#ifdef COIN_HAS_CPX

    const double *colsol=mip.getColSolution();
    nTotalNodes += CPXgetnodecnt(mip.getEnvironmentPtr(),mip.getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
    double new_lb = 0;
    int mipstat = CPXgetstat(mip.getEnvironmentPtr(),mip.getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
    //    if(mipstat ==
    if((mipstat != CPXMIP_OPTIMAL && mipstat != CPXMIP_OPTIMAL_TOL)) {
      int status = CPXgetbestobjval(mip.getEnvironmentPtr(),mip.getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL), &new_lb);
      lb=max(lb,new_lb);
      if(status)
        throw CoinError("Error in getting CPLEX best bound","IpCbcOACutGenerator2","siBestObj");
    }
    else
      lb = mip.getObjValue();

    CPXsetdblparam(mip.getEnvironmentPtr(),CPX_PARAM_CUTUP, mip.getInfinity());

#else
    //      if(!mip.isSolutionLimitReached())
    //          break;
    nTotalNodes += mip.getNodeCount();
    const double * colsol = mip.bestSolution();
    lb = mip.getBestPossibleObjValue();
    if(mip.bestSolution() == NULL)
      break;
#endif

    nMajorIt++;

    std::cout<<"Found new lower bound of : "<<lb<<std::endl;
    if ((ub - lb) < precision)
      break;
    bool same = 1;
    for(int i = 0 ; i < numIntCols ; i++) {
      double value;
      value = max(solver.getColLower()[inds[i]],colsol[inds[i]]);
      value = min(solver.getColUpper()[inds[i]],value);
      if(fabs(value - x[i])>0.0001) {
        same = 0;
      }
      x[i] = value;
      solver1.setColLower(inds[i], x[i]);
      solver1.setColUpper(inds[i], x[i]);
    }
    if(same) {
      std::cout<<"Converged on same solution"<<std::endl;
      lb += 1e10;
      break;
    }
    solver1.turnOnSolverOutput();
    solver1.initialSolve();
    if(solver1.isProvenOptimal()) {
      if (solution==NULL) solution = new double[numcols];
      CoinCopyN(solver1.getColSolution(), numcols, solution);
      std::cout<<pbName<<" OA found easible solution of value "
      <<solver1.getObjValue()<<" in "
      <<CoinCpuTime() -BeginTimeGLOB<<" seconds, "
      <<nMajorIt<<" major iterations, took"
      <<nTotalNodes<<" nodes."<<std::endl;

      OsiCuts cs;
      solver1.getOuterApproximation(cs,1, false, true);
      OsiSolverInterface::ApplyCutsReturnCode acRc;
      ub = min (solver1.getObjValue(), ub);
      acRc = solver.applyCuts(cs);
      // Print applyCuts return code
      std::cout <<cs.sizeCuts() <<" cuts were generated" <<std::endl;
      std::cout <<"  " <<acRc.getNumInconsistent() <<" were inconsistent" <<std::endl;
      std::cout <<"  " <<acRc.getNumInconsistentWrtIntegerModel()
      <<" were inconsistent for this problem" <<std::endl;
      std::cout <<"  " <<acRc.getNumInfeasible() <<" were infeasible" <<std::endl;
      std::cout <<"  " <<acRc.getNumIneffective() <<" were ineffective" <<std::endl;
      std::cout <<"  " <<acRc.getNumApplied() <<" were applied" <<std::endl;
      std::cout << std::endl << std::endl;

    }


    else {
      if (doFp)//Launch FP
      {
        //Reset NLP bounds
        for(int i = 0 ; i < numIntCols ; i++)
        {
          solver1.setColLower(inds[i], solver.getColLower()[inds[i]]);
          solver1.setColUpper(inds[i], solver.getColUpper()[inds[i]]);
        }
        ResolutionInformation infos;
        nTimesFPCalled++;
        int provenInfeas =0;

        //compute FP maxTime
        //	      firstIterationTime = 1.;
        double fpTime = min( params.maxTime_ - CoinCpuTime() + BeginTimeGLOB,120.);
        int fpMaxIter = 5;
        ub = min(ub, FP(solver1, solver,numIntCols, inds, x,fpTime, fpMaxIter,  infos, ub*(1-1e-03), provenInfeas, solution));
        if(provenInfeas == 1)
          lb = mip.getInfinity();
        if(provenInfeas == -1)
        {
          numNotFound++;
          // 		  if(numNotFound==10)
          // 		    doFp = 0;
        }
        // 	      if(infos.n_iterations >= 10)
        // 		doFp = 0;
        FP_infos+=infos;
      }
      else {
        std::cout<<"Adding feasibility cuts based on 1-norm of"
        <<"constraint satisfaction"<<std::endl;
        //	      solver1.getOuterApproximation(cs,1);
        solver1.getFeasibilityOuterApproximation( numIntCols, x, inds, cs, false, true);
        OsiSolverInterface::ApplyCutsReturnCode acRc;
        acRc = solver.applyCuts(cs);
        // Print applyCuts return code
        std::cout <<cs.sizeCuts() <<" cuts were generated" <<std::endl;
        std::cout <<"  " <<acRc.getNumInconsistent() <<" were inconsistent" <<std::endl;
        std::cout <<"  " <<acRc.getNumInconsistentWrtIntegerModel()
        <<" were inconsistent for this problem" <<std::endl;
        std::cout <<"  " <<acRc.getNumInfeasible() <<" were infeasible" <<std::endl;
        std::cout <<"  " <<acRc.getNumIneffective() <<" were ineffective" <<std::endl;
        std::cout <<"  " <<acRc.getNumApplied() <<" were applied" <<std::endl;
        std::cout << std::endl << std::endl;
      }
    }

  }
  std::string algoName=(doFp)?" enhanced OA":" OA";
  double time = CoinCpuTime() - BeginTimeGLOB;
  std::cout<<pbName<<algoName<<" found optimal solution of value "<<ub<<" (lower bound is "<<lb<<") in "<<time<<" seconds, "<<nMajorIt<<" major iterations, "<<std::endl;
  if(doFp)
  {
  std::cout<<nTimesInitFPCalled<<" calls to initial FP ( "
  <<nFPIterations <<" it, and "<< FPTime
  <<"sec)"
  <<std::endl
  <<nTimesFPCalled<<"subsequent calls to FP ( "
  <<FP_infos.n_iterations <<" it, and "<< FP_infos.time
  <<"sec)"<<std::endl;
  }
  std::cout<<"Nlp solve time : "<<nlpTime + FP_infos.nlp_time<<" B&B solve time : "<<bbTime + FP_infos.mip_time<<std::endl;
  delete[] inds;
  delete[] x;
  return 0;
}



/** Iterated feasibility pump.*/
int iteratedFP (AmplInterface& solver1, bool standAlone, 
                double * &solution)
{
  // Define a Solver which inherits from OsiClpsolverInterface -> OsiSolverInterface

  using namespace Ipopt;
  std::string pbName;
  solver1.getStrParam(OsiProbName,pbName);

#ifdef COIN_HAS_CPX

  OsiCpxSolverInterface solver;
  OsiCpxSolverInterface &mip = solver;
#else

  OsiClpSolverInterface solver;
#endif

  solver.messageHandler()->setLogLevel(0);

  //Setup timers since an nlp is solved to extract the linear relaxation
  double beginTime= CoinCpuTime();
  double bbTime = 0;
  double nlpTime = - beginTime;

  solver1.extractLinearRelaxation(solver, 1);
  nlpTime += CoinCpuTime();
  //Add extra variables for linearizing norm-1
  //get the number of 0-1 variables
  int numIntCols = 0;
  int numcols = solver1.getNumCols();
  int * inds = new int[numcols]; //indices of integer variables
  double * x = new double[numcols]; //to store values of integer variables later on

  for(int i = 0; i < numcols; i++) {
    if(solver1.isInteger(i)) {
      inds[numIntCols++] = i;
    }
  }
  //int k = 0;


  for(int i = 0;i<solver.getNumCols();i++)
    solver.setObjCoeff(i,0.);

  //Limits of the procedure
  //int nMaxNodes = 100000;

  //counters for iterations, nodes,...
  int nMajorIt = 0;
  int nTotalNodes = 0;


  double lb=-DBL_MAX;
  double ub=DBL_MAX;
  bool solved = 0;
  while(ub-lb>1e-05) {
    double ub = 1e100;
    bool feasible = 1;
    while(feasible && ( CoinCpuTime() - beginTime < params.maxTime_) ) {
      int nIt = 0;
      while(feasible & !solved && (CoinCpuTime() - beginTime < params.maxTime_)) {
        nIt++;
        for(int i = 0 ; i < numIntCols; i++) {
          solver.setObjCoeff(inds[i],1 - 2* solver1.getColSolution()[inds[i]]);
        }
        //change bound on alpha
        solver.setColUpper(solver.getNumCols()-1,ub);

#ifndef COIN_HAS_CPX
  
        solver.initialSolve();
        CbcModel mip(solver);
        CbcStrategyDefault defaultStrategy;
        mip.setStrategy(defaultStrategy);
        mip.solver()->messageHandler()->setLogLevel(0);

        //Add some heuristics to get feasible solutions
        mip.setMaximumSeconds(params.maxTime_ - CoinCpuTime() +beginTime);
        mip.setMaximumSolutions(1);
#else
        //	      CPXsetintparam(cpxSi->getEnvironmentPtr(),CPX_PARAM_NODELIM, nMaxNodes - nTotalNodes);
        CPXsetintparam(mip.getEnvironmentPtr(),CPX_PARAM_INTSOLLIM, 3);
        CPXsetdblparam(mip.getEnvironmentPtr(),CPX_PARAM_TILIM,
                       params.maxTime_ - CoinCpuTime() + beginTime);

#endif

        bbTime -= CoinCpuTime();
        mip.branchAndBound();
        bbTime += CoinCpuTime();
        //      mip.solver()->writeMps("oa");

        OsiCuts cs;
#ifdef COIN_HAS_CPX

        const double *colsol=mip.getColSolution();
        nTotalNodes += CPXgetnodecnt(mip.getEnvironmentPtr(),mip.getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
        int mipstat = CPXgetstat(mip.getEnvironmentPtr(),mip.getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
        if(mipstat == CPXMIP_INFEASIBLE || mipstat == CPXMIP_INForUNBD) {
          std::cout<<"Infeasible mip"<<std::endl;
          feasible = 0;
          break;
        }

#else
        if(mip.bestSolution() == NULL)
          break;
        if(mip.isNodeLimitReached()) {
          break;
        }
        //		if(!mip.isSolutionLimitReached())
        //			break;
        nTotalNodes += mip.getNodeCount();
        const double * colsol = mip.bestSolution();
#endif

        nMajorIt++;
        for(int i = 0 ; i < numIntCols ; i++) {
          x[i] = max(solver1.getColLower()[inds[i]],colsol[inds[i]]);
          x[i] = min(solver1.getColUpper()[inds[i]],x[i]);
          //				std::cout<<"Var "<<ind[i]<<" value "<<x[i]<<std::endl;
        }
        nlpTime -= CoinCpuTime();
        double dist = solver1.getFeasibilityOuterApproximation( numIntCols, x, inds, cs, false, true);
        nlpTime += CoinCpuTime();
        std::cout<<"Dist : "<<dist<<std::endl;
        if(dist < 1e-05)//do integer infeasibility check on variables
        {
          solved = 1;
          for(int i = 0 ; i < numcols && solved; i++)
          {
            if(solver1.isInteger(i)) {
              const double &value = solver1.getColSolution()[i];
              if(fabs( value- floor(value + 0.5)) > 1e-04) {
                solved = 0;
                std::cout<<"Variable "<<i<<" has integer infeasibility : "<<fabs( value- floor(value + 0.5))<<std::endl;
              }
            }
          }
          feasible = !solved;
        }

        OsiSolverInterface::ApplyCutsReturnCode acRc;
        acRc = solver.applyCuts(cs);
        // Print applyCuts return code
        std::cout <<cs.sizeCuts() <<" cuts were generated" <<std::endl;
        std::cout <<"  " <<acRc.getNumInconsistent() <<" were inconsistent" <<std::endl;
        std::cout <<"  " <<acRc.getNumInconsistentWrtIntegerModel()
        <<" were inconsistent for this problem" <<std::endl;
        std::cout <<"  " <<acRc.getNumInfeasible() <<" were infeasible" <<std::endl;
        std::cout <<"  " <<acRc.getNumIneffective() <<" were ineffective" <<std::endl;
        std::cout <<"  " <<acRc.getNumApplied() <<" were applied" <<std::endl;
        std::cout << std::endl << std::endl;
      }
      double time = CoinCpuTime() - beginTime;
      if(solved) {
        //Resolve the NLP with fixed variables and original objective function
        for(int i = 0; i < numIntCols; i++) {
          solver1.setColLower(inds[i], x[i]);
          solver1.setColUpper(inds[i], x[i]);
        }
        solver1.turnOnSolverOutput();
        solver1.initialSolve();
        if(solver1.isProvenOptimal()) {
          ub = min(solver1.getObjValue() * (1-1e-4), ub);
          if (solution==NULL) solution = new double[numcols];
          CoinCopyN(solver1.getColSolution(), numcols, solution);
          std::cout<<pbName<<" FP found easible solution of value "
          <<solver1.getObjValue()<<" in "
          <<CoinCpuTime() - beginTime<<" seconds, "
          <<nMajorIt<<" major iterations, took"
          <<nTotalNodes<<" nodes."<<std::endl;
          solved = 1;
          OsiCuts cs;
          solver1.getOuterApproximation(cs,1, NULL, true);
          OsiSolverInterface::ApplyCutsReturnCode acRc;
          acRc = solver.applyCuts(cs);
          for(int i = 0 ; i < numIntCols ; i++) {
            solver1.setColLower(inds[i], solver.getColLower()[inds[i]]);
            solver1.setColUpper(inds[i], solver.getColUpper()[inds[i]]);
          }
          //  solver1.initialSolve();
        }
        else
          std::cout<<pbName<<" FP converged in "<<time<<" seconds, "<<nMajorIt<<" major iterations, took"
          <<nTotalNodes<<" nodes, but solution seems infeasible"<<std::endl;

      }
      else if (feasible) {
        //restart from NLP-relaxation optimum
        solver1.initialSolve();
        nIt = 0;
      }
    }
      if(standAlone)
        break;
  }
  int returncode = 0;
  if(!standAlone)
  {
  //double time = CoinCpuTime() - beginTime;
  if(solved) {
    std::cout<<"iterated FP finished in : "<<CoinCpuTime() - beginTime
    <<", value of optimum "<<ub
    <<", "<<nMajorIt<<" major iterations, took"
    <<nTotalNodes<<" nodes."<<std::endl;
    returncode = 1;
  }
  else
    std::cout<<"iterated FP aborted on time limit elapsed time : "
    <<CoinCpuTime() - beginTime<<",  "<<nMajorIt
    <<" major iterations, took"
    <<nTotalNodes<<" nodes."<<std::endl;
    returncode = 0;
  }
  else returncode = 0;
  std::cout<<"Nlp solve time : "<<nlpTime
  <<" B&B solve time : "<<bbTime<<std::endl;
  delete[] inds;
  delete[] x;
  return returncode;
}



double FPGeneralIntegers(AmplInterface &nlp, OsiSolverInterface &linearModel,
                         int numIntCols, int * inds, double * vals, double maxTime,
                         ResolutionInformation& info, double ub, bool &provenInfeas)

{
  provenInfeas=0;
#ifdef COIN_HAS_CPX

  OsiCpxSolverInterface * cpxSi = dynamic_cast<OsiCpxSolverInterface *>
                                  (&linearModel);
  OsiCpxSolverInterface &mip = *cpxSi;
#endif

  int &nMajorIt = info.n_iterations;
  int &nTotalNodes = info.n_mip_nodes;
  int &nNlpIterations = info.n_nlp_iterations;
  double& bbTime = info.mip_time;
  double& nlpTime = info.nlp_time;
  double &time = info.time = -CoinCpuTime();
  double objValue = DBL_MAX;
  {
    OsiCuts cs;
    //First get the closest point to the current integer (NLP-infeasible) optimum
    nlpTime -= CoinCpuTime();
    double dist = nlp.getFeasibilityOuterApproximation( numIntCols, vals, inds, cs, false, true);
    nlpTime += CoinCpuTime();
    if(dist < 1e-06)//do integer infeasibility check on variables
    {
      //Something wrong shouldn't be
      std::cout<<"Feasibility subproblem has objective 0 while problem was claimed infeasible before"<<std::endl
      <<"I am confused, exiting with error"<<std::endl;
      throw -1;
    }

    OsiSolverInterface::ApplyCutsReturnCode acRc;

    //       linearModel.writeMps("test1");
    //       acRc = linearModel.applyCuts(cs);
    //       linearModel.writeMps("test2");

    // Print applyCuts return code
    std::cout <<cs.sizeCuts() <<" cuts were generated" <<std::endl;
    std::cout <<"  " <<acRc.getNumInconsistent() <<" were inconsistent" <<std::endl;
    std::cout <<"  " <<acRc.getNumInconsistentWrtIntegerModel()
    <<" were inconsistent for this problem" <<std::endl;
    std::cout <<"  " <<acRc.getNumInfeasible() <<" were infeasible" <<std::endl;
    std::cout <<"  " <<acRc.getNumIneffective() <<" were ineffective" <<std::endl;
    std::cout <<"  " <<acRc.getNumApplied() <<" were applied" <<std::endl;
    std::cout << std::endl << std::endl;
  }
  //Modify the linearModel for feasibility pump
  //save objective
  double * saveObj = new double[linearModel.getNumCols()];
  CoinCopyN(linearModel.getObjCoefficients(), linearModel.getNumCols(), saveObj);
  //Add extra variables for linearizing norm-1
  //get the number of 0-1 variables

  int origNumCols = linearModel.getNumCols();
  int origNumRows = linearModel.getNumRows();

  double * colLb = new double[2*numIntCols+1];
  double * colUb = new double[2*numIntCols+1];
  double * obj = new double[2*numIntCols];
  int k = 0;
  CoinPackedVectorBase ** emptyCols = new CoinPackedVectorBase*[2*numIntCols+1];
  for(int i = 0; i < nlp.getNumCols(); i++) {
    if(nlp.isInteger(i)) {
      colLb[k] = colLb[k+1] = 0.;
      colUb[k] = colUb[k+1] = nlp.getColUpper()[i] - nlp.getColLower()[i];
      obj[k] = obj[k+1] = 1.;
      emptyCols[k] = new CoinPackedVector;
      emptyCols[k + 1] = new CoinPackedVector;
      k+=2;
    }
  }

  linearModel.addCols(2*numIntCols, emptyCols, colLb, colUb, obj);
  k = 0;
  //add the constraint x_i - x^+_i + x^-_i = 0
  for(int i = 0; i < nlp.getNumCols(); i++) {
    if(nlp.isInteger(i)) {
      CoinPackedVector * rowToAdd = (CoinPackedVector *) emptyCols[k];
      rowToAdd->insert(i,1.);
      rowToAdd->insert(origNumCols + 2*k ,-1.);
      rowToAdd->insert(origNumCols + 2*k + 1,1.);
      colUb[k] = 0.;
      k++;
    }
  }
  //int addCutOff = 0;
  //if ub is given add a "cutoff" row
  if(ub < 1e100) {
    //change bound on alpha
    linearModel.setColUpper(origNumCols-1,ub);
    //     addCutOff=1;
    //     CoinPackedVector * rowToAdd = (CoinPackedVector *) emptyCols[k];
    //     rowToAdd->insert(origNumCols-1,1.);
    //     colUb[k] = ub;
    //     colLb[k] = -linearModel.getInfinity();
  }
  for(int i = 0;i<origNumCols;i++)
    linearModel.setObjCoeff(i,0.);
  linearModel.addRows(numIntCols, emptyCols, colLb, colUb);
  //done
  int nAddedCuts = 0;

  bool solved = 0;
  bool numericFailure = 0;
  int numScaleIncrease = 0;
  while(!solved && (time + CoinCpuTime() < maxTime) && !numericFailure && numScaleIncrease < 5 ) {
    for(int i = 0; i < numIntCols; i++) {
      linearModel.setRowBounds(origNumRows + i, nlp.getColSolution()[inds[i]], nlp.getColSolution()[inds[i]]);
    }
#ifdef COIN_HAS_CPX
    //	      CPXsetintparam(cpxSi->getEnvironmentPtr(),CPX_PARAM_NODELIM, nMaxNodes - nTotalNodes);
    CPXsetintparam(cpxSi->getEnvironmentPtr(),CPX_PARAM_INTSOLLIM, 1);
    CPXsetdblparam(cpxSi->getEnvironmentPtr(),CPX_PARAM_TILIM, maxTime - CoinCpuTime() + time);

#else

    linearModel.initialSolve();
    CbcModel mip(linearModel);
    //        mip.writeMps("oa");
    CbcStrategyDefault defaultStrategy;
    mip.setStrategy(defaultStrategy);
    mip.solver()->messageHandler()->setLogLevel(0);

    mip.setMaximumSeconds(maxTime - CoinCpuTime() + time);
    mip.setMaximumSolutions(1);
#endif

    bbTime -= CoinCpuTime();
    mip.branchAndBound();
    bbTime += CoinCpuTime();
    //      mip.writeMps("oa1");
    OsiCuts cs;
    nMajorIt++;

#ifdef COIN_HAS_CPX

    const double * colsol = mip.getColSolution();
    //int nNodes = 0;
    nTotalNodes += CPXgetnodecnt(cpxSi->getEnvironmentPtr(),cpxSi->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
    int mipstat = CPXgetstat(mip.getEnvironmentPtr(),mip.getLpPtr(OsiCpxSolverInterface::KEEPCACHED_ALL));
    if(mipstat == CPXMIP_INFEASIBLE) {
      provenInfeas=1;
      return 1e75;
    }
#else
    if(mip.bestSolution() == NULL)
      break;
    if(mip.isNodeLimitReached()) {
      break;
    }
    //      if(!mip.isSolutionLimitReached())
    //          break;
    nTotalNodes += mip.getNodeCount();
    const double * colsol = mip.bestSolution();
#endif

    for(int i = 0 ; i < numIntCols ; i++) {
      vals[i] = floor(colsol[inds[i]] + 0.5);
      vals[i] = max(mip.getColLower()[inds[i]],vals[i]);
      vals[i] = min(mip.getColUpper()[inds[i]],vals[i]);
    }
    nlpTime -= CoinCpuTime();
    double dist = nlp.getFeasibilityOuterApproximation( numIntCols, vals, inds, cs, false, true);
    nlpTime += CoinCpuTime();
    nNlpIterations += nlp.getIterationCount();
    if(!nlp.isProvenOptimal()) {
      std::cout<<"?????"<<std::endl;
      throw -1;
    }
    nlpTime += CoinCpuTime();
    if(dist < 1e-04)//do integer infeasibility check on variables
    {
      solved = 1;
      double norm_inf = DBL_MAX;
      for(int i = 0 ; i < nlp.getNumCols() && solved; i++)
      {
        if(nlp.isInteger(i)) {
          const double &value = nlp.getColSolution()[i];
          double IIf = fabs( value- floor(value + 0.5));
          norm_inf = min(norm_inf, IIf);
          if(fabs( value- floor(value + 0.5)) > 1e-04) {
            std::cout<<"Variable "<<i<<" has integer infeasibility : "<<IIf<<std::endl;
            solved=0;
            numericFailure = 1;
#if 0

            numScaleIncrease++;
            nlp.fpnlp()->setObjectiveScaling(10* nlp.fpnlp()->getObjectiveScaling());
            for(int i = 0 ; i < numIntCols ; i++) {
              nlp.setColLower(inds[i], linearModel.getColLower()[inds[i]]);
              nlp.setColUpper(inds[i], linearModel.getColUpper()[inds[i]]);
            }
#endif

          }
        }
      }
      std::cout<<"Found a solution with maximal integer infeasibility of "<<norm_inf<<std::endl;
    }

    OsiSolverInterface::ApplyCutsReturnCode acRc;
    std::cout<<"Cut violation :"<<cs.rowCut(0).violated(colsol)
    <<std::endl;
    //       cs.printCuts();

    linearModel.writeMps("test1");
    acRc = linearModel.applyCuts(cs);
    linearModel.writeMps("test2");
    nAddedCuts += cs.sizeRowCuts();
    // Print applyCuts return code
    std::cout <<cs.sizeCuts() <<" cuts were generated" <<std::endl;
    std::cout <<"  " <<acRc.getNumInconsistent() <<" were inconsistent" <<std::endl;
    std::cout <<"  " <<acRc.getNumInconsistentWrtIntegerModel()
    <<" were inconsistent for this problem" <<std::endl;
    std::cout <<"  " <<acRc.getNumInfeasible() <<" were infeasible" <<std::endl;
    std::cout <<"  " <<acRc.getNumIneffective() <<" were ineffective" <<std::endl;
    std::cout <<"  " <<acRc.getNumApplied() <<" were applied" <<std::endl;
    std::cout << std::endl << std::endl;


    if(solved) {
      //Set warm start point to the last point found (which is feasible for this relaxation)
      nlp.setColSolution(nlp.getColSolution());
      nlp.setRowPrice(nlp.getRowPrice());
      nlp.solver()->enableWarmStart();
      //Resolve the NLP with fixed variables and original objective function
      for(int i = 0; i < numIntCols; i++) {
        nlp.setColLower(inds[i], vals[i]);
        nlp.setColUpper(inds[i], vals[i]);
      }
      //nlp.turnOnSolverOutput();
      nlp.initialSolve();
      if(nlp.isProvenOptimal()) {
        OsiCuts cs;
        nlp.getOuterApproximation(cs,1, NULL, true);
        linearModel.applyCuts(cs);
        nAddedCuts += cs.sizeRowCuts();
        std::cout<<" FP found easible solution of value "<<nlp.getObjValue()<<" in "<<time + CoinCpuTime()<<" seconds, "<<nMajorIt<<" major iterations, took"
        <<nTotalNodes<<" nodes."<<std::endl;
        objValue = nlp.getObjValue();
      }
      else {
        std::cout<<" FP converged in "<<time<<" seconds, "<<nMajorIt<<" major iterations, took"
        <<nTotalNodes<<" nodes, but solution seems infeasible"<<std::endl;
        std::cout<<"Increasing scaling of objective and restarting"<<std::endl;
        solved = 0;
        numScaleIncrease++;
        numericFailure = 1;
#if 0

        nlp.fpnlp()->setObjectiveScaling(10* nlp.fpnlp()->getObjectiveScaling());
        for(int i = 0 ; i < numIntCols ; i++) {
          nlp.setColLower(inds[i], linearModel.getColLower()[inds[i]]);
          nlp.setColUpper(inds[i], linearModel.getColUpper()[inds[i]]);
        }
#endif

      }

    }

  }
  time+=CoinCpuTime();

  if(!solved) {
    if(numericFailure)
      std::cout<<"FP aborted because of a numberical failure : "<<time<<",  "<<nMajorIt<<" major iterations, took"
      <<nTotalNodes<<" nodes."<<std::endl;
    else
      std::cout<<"FP aborted on time limit elapsed time : "<<time<<",  "<<nMajorIt<<" major iterations, took"
      <<nTotalNodes<<" nodes."<<std::endl;

  }
  std::cout<<"Nlp solve time : "<<nlpTime<<" B&B solve time : "<<bbTime<<std::endl;


  //erase the extra columns and extra rows added to linearModel and resore the objective function
  int * colToDelete = new int[2*numIntCols];
  int * rowToDelete = new int[numIntCols];
  for(int i = 0 ; i < numIntCols ; i++) {
    colToDelete[2*i] = origNumCols + 2*i;
    colToDelete[2*i +1] = origNumCols + 2*i + 1;
    rowToDelete[i] = origNumRows + i;
  }
  //   if(addCutOff)
  //     rowToDelete[numIntCols] = linearModel.getNumCols();
  linearModel.deleteCols(2*numIntCols, colToDelete);
  linearModel.deleteRows(numIntCols, rowToDelete);
  //Check that size are correct
  //  mip.writeMps("oa2");
  assert(linearModel.getNumCols()==origNumCols);
  assert(linearModel.getNumRows()==origNumRows + nAddedCuts);
  //    restore objective
  for(int i = 0 ; i < origNumCols; i++)
    linearModel.setObjCoeff(i, saveObj[i]);
  delete [] saveObj;
  return objValue;
}

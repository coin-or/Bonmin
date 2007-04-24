// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/13/2007

#include "OsiClpSolverInterface.hpp"

#include "BonBonminSetup.hpp"
#include "BonCurvBranching.hpp"
#include "BonQPStrongBranching.hpp"
#include "BonLpStrongBranching.hpp"

//OA machinery
#include "BonDummyHeuristic.hpp"
#include "BonOACutGenerator2.hpp"
#include "BonOaFeasChecker.hpp"
#include "BonOaNlpOptim.hpp"
#include "BonEcpCuts.hpp"

#include "BonCbcNode.hpp"



//MILP cuts
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglTwomir.hpp"
#include "CglPreProcess.hpp"
#include "CglLandP.hpp"
#include "CglRedSplit.hpp"

namespace Bonmin{
  BonminSetup::BonminSetup():BabSetupBase(){
  }
  
  BonminSetup::BonminSetup(BasicSetup& b, Ipopt::SmartPtr<TMINLP> tminlp):
  BabSetupBase(b, tminlp){
  }
  
  BonminSetup::BonminSetup(const OsiTMINLPInterface& nlp):
  BabSetupBase(nlp){
  }
  
  BonminSetup::BonminSetup(const BonminSetup &other):BabSetupBase(other){
  }
  
  void BonminSetup::registerAllOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions){
    BabSetupBase::registerAllOptions(roptions);
    /* Branching options.*/
    LpStrongBranching::registerOptions(roptions);
    
    /* Outer Approximation options.*/
    OACutGenerator2::registerOptions(roptions);
    EcpCuts::registerOptions(roptions);
    OaNlpOptim::registerOptions(roptions);
    
    BonCbcFullNodeInfo::registerOptions(roptions);
    
    registerMilpCutGenerators(roptions);
    
    
    roptions->SetRegisteringCategory("bonmin options: Branching strategies options");
    roptions->AddStringOption2("sos_constraints",
                               "Wether or not to activate SOS constraints.",
                               "enable",
                               "enable","",
                               "disable","",
                               "(only type 1 SOS are supported at the moment)");
    
    roptions->AddStringOption9("varselect_stra",
                               "Chooses variable selection strategy",
                               "strong-branching",
                               "most-fractional", "Choose most fractional variable",
                               "strong-branching", "Perform strong branching",
                               "reliability-branching", "Use reliability branching",
                               "curvature-estimator", "Use curvature estimation to select branching variable",
                               "qp-strong-branching", "Perform strong branching with QP approximation",
                               "lp-strong-branching", "Perform strong branching with LP approximation",
                               "nlp-strong-branching", "Perform strong branching with NLP approximation",
                               "osi-simple", "Osi method to do simple branching",
                               "osi-strong", "Osi method to do strong branching","");
    
    
  }
  /** Register all the Bonmin options.*/
  void 
  BonminSetup::registerOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions){
    registerAllOptions(roptions);
  }
    
  /** Initialize, read options and create appropriate bonmin setup using initialized tminlp.*/
  void 
  BonminSetup::initializeBonmin(Ipopt::SmartPtr<TMINLP> tminlp){
     defaultBasicOptions();

    initialize(tminlp);
    int ival;
    options_->GetEnumValue("algorithm", ival,"bonmin.");
    BaseOptions::Algorithm algo = (BaseOptions::Algorithm) ival;
    if(algo == BaseOptions::B_BB)
      initializeBBB();
    else
      initializeBHyb();}
  

    void
    BonminSetup::defaultBasicOptions(){
      if(GetRawPtr(options_) != NULL && GetRawPtr(roptions_) != NULL &&  GetRawPtr(journalist_) != NULL) return;
      BasicSetup b;
      Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions = b.roptions();
      BonminSetup::registerAllOptions(roptions);
      b.Initialize("bonmin.opt");
      options_ = b.options();
      roptions_ = b.roptions();
      journalist_ = b.journalist();
    }
    
    
  
  /** Initialize, read options and create appropriate bonmin setup using initialized tminlp.*/
  void 
  BonminSetup::initializeBonmin(const OsiTMINLPInterface &nlpSi){
    defaultBasicOptions();
    initialize(nlpSi);
    int ival;
    options_->GetEnumValue("algorithm", ival,"bonmin.");
    BaseOptions::Algorithm algo = (BaseOptions::Algorithm) ival;
    if(algo == BaseOptions::B_BB)
      initializeBBB();
    else
      initializeBHyb();
  }
  
  /** Register standard MILP cut generators. */
  void 
  BonminSetup::registerMilpCutGenerators(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions){
    roptions->SetRegisteringCategory("bonmin options for MILP cutting planes");
    
    roptions->AddLowerBoundedIntegerOption("Gomory_cuts",
                                           "Frequency k (in terms of nodes) for generating Gomory cuts in branch-and-cut.",
                                           -100,-5,
                                           "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
                                           "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
                                           "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
    roptions->AddLowerBoundedIntegerOption("probing_cuts",
                                           "Frequency (in terms of nodes) for generating probing cuts in branch-and-cut",
                                           -100,-5,
                                           "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
                                           "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
                                           "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
    
    roptions->AddLowerBoundedIntegerOption("cover_cuts",
                                           "Frequency (in terms of nodes) for generating cover cuts in branch-and-cut",
                                           -100,-5,
                                           "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
                                           "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
                                           "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
    
    roptions->AddLowerBoundedIntegerOption("mir_cuts",
                                           "Frequency (in terms of nodes) for generating MIR cuts in branch-and-cut",
                                           -100,-5,
                                           "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
                                           "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
                                           "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
    roptions->AddLowerBoundedIntegerOption("2mir_cuts",
                                           "Frequency (in terms of nodes) for generating 2-MIR cuts in branch-and-cut",
                                           -100,0,
                                           "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
                                           "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
                                           "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
    roptions->AddLowerBoundedIntegerOption("flow_covers_cuts",
                                           "Frequency (in terms of nodes) for generating flow cover cuts in branch-and-cut",
                                           -100,-5,
                                           "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
                                           "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
                                           "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
    roptions->AddLowerBoundedIntegerOption("lift_and_project_cuts",
                                           "Frequency (in terms of nodes) for generating lift-and-project cuts in branch-and-cut",
                                           -100,0,
                                           "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
                                           "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
                                           "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
    roptions->AddLowerBoundedIntegerOption("reduce_and_split_cuts",
                                           "Frequency (in terms of nodes) for generating reduce-and-split cuts in branch-and-cut",
                                           -100,0,
                                           "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
                                           "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
                                           "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
    roptions->AddLowerBoundedIntegerOption("clique_cuts",
                                           "Frequency (in terms of nodes) for generating clique cuts in branch-and-cut",
                                           -100,-5,
                                           "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
                                           "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
                                           "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
  }
  /** Add milp cut generators according to options.*/
  void 
  BonminSetup::addMilpCutGenerators(){
    int freq;
    options_->GetIntegerValue("Gomory_cuts", freq,"bonmin.");
    if(freq){
      CuttingMethod cg;
      cg.frequency = freq;
      CglGomory * gom = new CglGomory;
      cg.cgl = gom;
      gom->setLimitAtRoot(512);
      gom->setLimit(50);
      cg.id = "Mixed Integer Gomory";
      cutGenerators_.push_back(cg);
    }
    options_->GetIntegerValue("probing_cuts",freq,"bonmin.");
    if(freq){
      CuttingMethod cg;
      cg.frequency = freq;
      CglProbing * probe = new CglProbing;
      cg.cgl = probe;
      probe->setUsingObjective(1);
      probe->setMaxPass(3);
      probe->setMaxPassRoot(3);
      // Number of unsatisfied variables to look at
      probe->setMaxProbe(10);
      probe->setMaxProbeRoot(50);
      // How far to follow the consequences
      probe->setMaxLook(10);
      probe->setMaxLookRoot(50);
      probe->setMaxLookRoot(10);
      // Only look at rows with fewer than this number of elements
      probe->setMaxElements(200);
      probe->setRowCuts(3);
      cg.id = "Probing";
      cutGenerators_.push_back(cg);
    }
    options_->GetIntegerValue("mir_cuts",freq,"bonmin.");
    if(freq){
      CuttingMethod cg;
      cg.frequency = freq;
      CglMixedIntegerRounding2 * mir = new CglMixedIntegerRounding2;
      cg.cgl = mir;
      cg.id = "Mixed Integer Rounding";
      cutGenerators_.push_back(cg);
      
      
    }
    options_->GetIntegerValue("2mir_cuts",freq,"bonmin.");
    if(freq){
      CuttingMethod cg;
      cg.frequency = freq;
      CglTwomir * mir2 = new CglTwomir;
      cg.cgl = mir2;
      cg.id = "2-MIR";
      cutGenerators_.push_back(cg);
    }      
    options_->GetIntegerValue("cover_cuts",freq,"bonmin.");
    if(freq){
      CuttingMethod cg;
      cg.frequency = freq;
      CglKnapsackCover * cover = new CglKnapsackCover;
      cg.cgl = cover;
      cg.id = "Cover";
      cutGenerators_.push_back(cg);
    }
    
    options_->GetIntegerValue("clique_cuts",freq,"bonmin.");
    if(freq){
      CuttingMethod cg;
      cg.frequency = freq;
      CglClique * clique = new CglClique;
      clique->setStarCliqueReport(false);
      clique->setRowCliqueReport(false);
      clique->setMinViolation(0.1);
      
      cg.cgl = clique;
      cg.id = "Clique";
      cutGenerators_.push_back(cg);      
    }
    options_->GetIntegerValue("flow_covers_cuts",freq,"bonmin.");
    if(freq){
      CuttingMethod cg;
      cg.frequency = freq;
      CglFlowCover * flow = new CglFlowCover;
      cg.cgl = flow;
      cg.id = "Flow Covers";
      cutGenerators_.push_back(cg);
    }
    options_->GetIntegerValue("lift_and_project_cuts",freq,"bonmin.");
    if(freq){
      CuttingMethod cg;
      cg.frequency = freq;
      CglLandP * landp = new CglLandP;
      cg.cgl = landp;
      cg.id = "Lift-and-Project";
      cutGenerators_.push_back(cg);
    }
    options_->GetIntegerValue("reduce_and_split_cuts",freq,"bonmin.");
    if(freq){
      CuttingMethod cg;
      cg.frequency = freq;
      CglRedSplit * rands = new CglRedSplit;
      cg.cgl = rands;
      cg.id = "Reduce-and-Split";
      cutGenerators_.push_back(cg);
    }
  }
  
  
  void 
  BonminSetup::initializeBBB(){
    linearSolver_ = nonlinearSolver_;
    nonlinearSolver_->ignoreFailures();
    OsiBabSolver extraStuff(2);
    linearSolver_->setAuxiliaryInfo(&extraStuff);
    
    intParam_[BabSetupBase::SpecialOption] = 16;
    intParam_[BabSetupBase::MinReliability] = 0;
    if(!options_->GetIntegerValue("number_strong_branch",intParam_[NumberStrong],"bonmin.")){
      intParam_[BabSetupBase::NumberStrong] = 0;
    }
    int varSelection;
    options_->GetEnumValue("varselect_stra",varSelection,"bonmin.");
    
    if(varSelection == OsiTMINLPInterface::CURVATURE_ESTIMATOR){
      linearSolver_->findIntegersAndSOS(false);
      BonCurvBranching * chooseVariable = new BonCurvBranching(nonlinearSolver_);
      branchingMethod_ = chooseVariable;      
    }
    else if(varSelection == OsiTMINLPInterface::QP_STRONG_BRANCHING){
      linearSolver_->findIntegersAndSOS(false);
      BonQPStrongBranching*  chooseVariable = new BonQPStrongBranching(nonlinearSolver_);
      branchingMethod_ = chooseVariable;
    }
    else if(varSelection == OsiTMINLPInterface::LP_STRONG_BRANCHING){
      linearSolver_->findIntegersAndSOS(false);
      LpStrongBranching * choose = new LpStrongBranching(nonlinearSolver_);
      choose->setMaxCuttingPlaneIter(intParam_[BabSetupBase::NumEcpRoundsStrong]);
      branchingMethod_ = choose;
    }
    else if(varSelection == OsiTMINLPInterface::NLP_STRONG_BRANCHING){
      const bool solve_nlp = true;
      linearSolver_->findIntegersAndSOS(false);
      BonQPStrongBranching * choose =  new BonQPStrongBranching(nonlinearSolver_, solve_nlp);
      branchingMethod_ = choose;
    }
    else if(varSelection == OsiTMINLPInterface::OSI_SIMPLE){
      linearSolver_->findIntegersAndSOS(false);
      branchingMethod_ = new OsiChooseVariable(nonlinearSolver_);
    }
    else if(varSelection == OsiTMINLPInterface::OSI_STRONG){
      linearSolver_->findIntegersAndSOS(false);
      branchingMethod_ = new OsiChooseStrong(nonlinearSolver_);
    }
  }  
  
  void 
  BonminSetup::initializeBHyb()
{
    /* Create linear solver */
    linearSolver_ = new OsiClpSolverInterface;
    int lpLogLevel;
    options_->GetIntegerValue("lp_log_level",lpLogLevel,"bonmin.");
    linearSolver_->messageHandler()->setLogLevel(lpLogLevel);
    nonlinearSolver_->extractLinearRelaxation(*linearSolver_);
    
    int ival;
    options_->GetEnumValue("algorithm", ival,"bonmin.");
    BaseOptions::Algorithm algo = (BaseOptions::Algorithm) ival;
    if(algo == BaseOptions::B_OA){
      options_->SetNumericValue("oa_dec_time_limit",DBL_MAX, true, true);
      options_->SetNumericValue("nlp_solve_frequency", 0, true, true);
      intParam_[BabLogLevel] = 0;
    }
    else if (algo==BaseOptions::B_QG) {
      options_->SetNumericValue("oa_dec_time_limit",0, true, true);
      options_->SetNumericValue("nlp_solve_frequency", 0, true, true);
    }
    //#define GREAT_STUFF_FOR_ANDREAS
#ifdef GREAT_STUFF_FOR_ANDREAS
    printf("ToDo: Clean me up in Bab::branchAndBound\n");
    OsiCuts cuts;
    nonlinearSolver_->getOuterApproximation(cuts, true, NULL, true);
    si->applyCuts(cuts);
#endif
    // say bound dubious, does cuts at solution
    OsiBabSolver * extraStuff = new OsiBabSolver(3);
    linearSolver_->setAuxiliaryInfo(extraStuff);
    delete extraStuff;
    
    int varSelection;
    options_->GetEnumValue("varselect_stra",varSelection,"bonmin.");
    if(varSelection > OsiTMINLPInterface::RELIABILITY_BRANCHING){
      std::cout<<"Variable selection stragey not available with oa branch-and-cut."<<std::endl;
    }
    /* Populate cut generation and heuristic procedures.*/
    options_->GetIntegerValue("nlp_solve_frequency",ival,"bonmin.");
    if(ival != 0){
      CuttingMethod cg;
      cg.frequency = ival;
      OaNlpOptim * nlpsol = new OaNlpOptim(*this);
      cg.cgl = nlpsol;
      cg.id="NLP solution based oa cuts";
      cutGenerators_.push_back(cg);
    }
    
    options_->GetIntegerValue("filmint_ecp_cuts",ival, "bonmin.");
    if(ival != 0){
      CuttingMethod cg;
      cg.frequency = ival;
      EcpCuts * ecp = new EcpCuts(*this);
      cg.cgl = ecp;
      cg.id = "Ecp cuts";
      cutGenerators_.push_back(cg);
    }
    
    if (algo!=BaseOptions::B_QG)
      addMilpCutGenerators();
    
    double oaTime;
    options_->GetNumericValue("oa_dec_time_limit",oaTime,"bonmin.");
    if(oaTime > 0.)
    {
      CuttingMethod cg;
      cg.frequency = ival;
      OACutGenerator2 * oa = new OACutGenerator2(*this);
      cg.cgl = oa;
      cg.id = "Outer Approximation decomposition.";
      cutGenerators_.push_back(cg);
      
    }
    
    {
      CuttingMethod cg;
      cg.frequency = 1;
      OaFeasibilityChecker * oa = new OaFeasibilityChecker(*this);
      cg.cgl = oa;
      cg.id = "Outer Approximation feasibility check.";
      cg.atSolution = 1;
      cutGenerators_.push_back(cg);
    }
}
}/* end namespace Bonmin*/


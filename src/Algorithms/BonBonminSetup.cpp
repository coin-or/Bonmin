// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/13/2007

#include "BonminConfig.h"
#include "OsiClpSolverInterface.hpp"

#include "BonBonminSetup.hpp"
#ifdef BONMIN_CURVATURE_BRANCHING
#include "BonCurvBranchingSolver.hpp"
#endif
#include "BonChooseVariable.hpp"
#include "BonRandomChoice.hpp"
#include "BonDiver.hpp"
#include "BonQpBranchingSolver.hpp"
#include "BonLpBranchingSolver.hpp"

//OA machinery
#include "BonDummyHeuristic.hpp"
#include "BonOACutGenerator2.hpp"
#include "BonFpForMinlp.hpp"
#include "BonOaFeasChecker.hpp"
#include "BonOaNlpOptim.hpp"
#include "BonEcpCuts.hpp"

#include "BonCbcNode.hpp"
#ifdef COIN_HAS_FILTERSQP
# include "BonFilterSolver.hpp"
#endif

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
#include "BonLinearCutsGenerator.hpp"

#include "BonFixAndSolveHeuristic.hpp"
#include "BonDummyPump.hpp"
#include "BonPumpForMinlp.hpp"
#include "BonHeuristicRINS.hpp"
#include "BonHeuristicLocalBranching.hpp"
#include "BonHeuristicFPump.hpp"
#include "BonHeuristicDiveFractional.hpp"
#include "BonHeuristicDiveVectorLength.hpp"
#include "BonHeuristicDiveMIPFractional.hpp"
#include "BonHeuristicDiveMIPVectorLength.hpp"
#include "BonMilpRounding.hpp"
//#include "BonInnerApproximation.hpp"
namespace Bonmin
{
  BonminSetup::BonminSetup(const CoinMessageHandler * handler):BabSetupBase(handler),algo_(Dummy)
  {}

  BonminSetup::BonminSetup(const BonminSetup &other):BabSetupBase(other),
      algo_(other.algo_)
  {}

  BonminSetup::BonminSetup(const BonminSetup &other,
                           OsiTMINLPInterface &nlp):
      BabSetupBase(other, nlp),
      algo_(other.algo_)
  {
    if(algo_ != B_BB){
      assert(continuousSolver_ == NULL);
      continuousSolver_ = new OsiClpSolverInterface;
      int lpLogLevel;
      options_->GetIntegerValue("lp_log_level",lpLogLevel,prefix_.c_str());
      if(messageHandler_)
        continuousSolver_->passInMessageHandler(messageHandler_);
      continuousSolver_->messageHandler()->setLogLevel(lpLogLevel);

      nonlinearSolver_->extractLinearRelaxation(*continuousSolver_);
      // say bound dubious, does cuts at solution
      OsiBabSolver * extraStuff = new OsiBabSolver(3);
      continuousSolver_->setAuxiliaryInfo(extraStuff);
      delete extraStuff;
    }
  }
  BonminSetup::BonminSetup(const BonminSetup &other,
                           OsiTMINLPInterface &nlp,
                           const std::string &prefix):
    BabSetupBase(other, nlp, prefix),
    algo_(Dummy)
  {
   algo_ = getAlgorithm();
    if (algo_ == B_BB)
      initializeBBB();
    else
      initializeBHyb(true);
  }
  void BonminSetup::registerAllOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    BabSetupBase::registerAllOptions(roptions);

    /* Outer Approximation options.*/
    OACutGenerator2::registerOptions(roptions);
    OaFeasibilityChecker::registerOptions(roptions);
    MinlpFeasPump::registerOptions(roptions);
    EcpCuts::registerOptions(roptions);
    OaNlpOptim::registerOptions(roptions);
    SubMipSolver::registerOptions(roptions);


    BonCbcFullNodeInfo::registerOptions(roptions);


    registerMilpCutGenerators(roptions);


    /** Heursitics.*/
    LocalSolverBasedHeuristic::registerOptions(roptions);
    FixAndSolveHeuristic::registerOptions(roptions);
    DummyPump::registerOptions(roptions);
    MilpRounding::registerOptions(roptions);
    PumpForMinlp::registerOptions(roptions);
    HeuristicRINS::registerOptions(roptions);
    HeuristicLocalBranching::registerOptions(roptions);
    HeuristicFPump::registerOptions(roptions);
    HeuristicDiveFractional::registerOptions(roptions);
    HeuristicDiveVectorLength::registerOptions(roptions);
    HeuristicDiveMIPFractional::registerOptions(roptions);
    HeuristicDiveMIPVectorLength::registerOptions(roptions);

    roptions->SetRegisteringCategory("Algorithm choice", RegisteredOptions::BonminCategory);
    roptions->AddStringOption6("algorithm",
        "Choice of the algorithm.",
        "B-BB",
        "B-BB","simple branch-and-bound algorithm,",
        "B-OA","OA Decomposition algorithm,",
        "B-QG","Quesada and Grossmann branch-and-cut algorithm,",
        "B-Hyb","hybrid outer approximation based branch-and-cut,",
        "B-Ecp","ECP cuts based branch-and-cut a la FilMINT.",
        "B-iFP","Iterated Feasibility Pump for MINLP.",
        "This will preset some of the options of bonmin depending on the algorithm choice."
                              );
    roptions->setOptionExtraInfo("algorithm",127);


  }

  /** Register all the Bonmin options.*/
  void
  BonminSetup::registerOptions()
  {
    registerAllOptions(roptions_);
  }

  /** Initialize, read options and create appropriate bonmin setup using initialized tminlp.*/
  void
  BonminSetup::initialize(Ipopt::SmartPtr<TMINLP> tminlp, bool createContinuousSolver /*= false*/)
  {

    use(tminlp);
    BabSetupBase::gatherParametersValues(options_);
    algo_ = getAlgorithm();
    if (algo_ == B_BB)
      initializeBBB();
    else
      initializeBHyb(createContinuousSolver);
  }

  /** Initialize, read options and create appropriate bonmin setup using initialized tminlp.*/
  void
  BonminSetup::initialize(const OsiTMINLPInterface &nlpSi, bool createContinuousSolver /*= false*/)
  {
    use(nlpSi);
    BabSetupBase::gatherParametersValues(options_);
    Algorithm algo = getAlgorithm();
    if (algo == B_BB)
      initializeBBB();
    else
      initializeBHyb(createContinuousSolver);
  }

  /** Register standard MILP cut generators. */
  void
  BonminSetup::registerMilpCutGenerators(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("MILP cutting planes in hybrid algorithm", RegisteredOptions::BonminCategory);

    roptions->AddLowerBoundedIntegerOption("Gomory_cuts",
        "Frequency (in terms of nodes) for generating Gomory cuts in branch-and-cut.",
        -100,-5,
        "If $k > 0$, cuts are generated every $k$ nodes, if $-99 < k < 0$ cuts are generated every $-k$ nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if $k=-99$ generate cuts only at the root node, if $k=0$ or $100$ do not generate cuts.");
    roptions->setOptionExtraInfo("Gomory_cuts",119);
#if 0
    roptions->AddBoundedIntegerOption("probing_cuts",
        "Frequency (in terms of nodes) for generating probing cuts in branch-and-cut",
        0,0,0,
        "If $k > 0$, cuts are generated every $k$ nodes, if $-99 < k < 0$ cuts are generated every $-k$ nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if $k=-99$ generate cuts only at the root node, if $k=0$ or $100$ do not generate cuts.");
    roptions->setOptionExtraInfo("probing_cuts",0);
#endif
    roptions->AddLowerBoundedIntegerOption("cover_cuts",
        "Frequency (in terms of nodes) for generating cover cuts in branch-and-cut",
        -100,0,
        "If $k > 0$, cuts are generated every $k$ nodes, if $-99 < k < 0$ cuts are generated every $-k$ nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if $k=-99$ generate cuts only at the root node, if $k=0$ or $100$ do not generate cuts.");
    roptions->setOptionExtraInfo("cover_cuts",119);

    roptions->AddLowerBoundedIntegerOption("mir_cuts",
        "Frequency (in terms of nodes) for generating MIR cuts in branch-and-cut",
        -100,-5,
        "If $k > 0$, cuts are generated every $k$ nodes, if $-99 < k < 0$ cuts are generated every $-k$ nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if $k=-99$ generate cuts only at the root node, if $k=0$ or $100$ do not generate cuts.");
    roptions->setOptionExtraInfo("mir_cuts",119);
    roptions->AddLowerBoundedIntegerOption("2mir_cuts",
        "Frequency (in terms of nodes) for generating 2-MIR cuts in branch-and-cut",
        -100,0,
        "If $k > 0$, cuts are generated every $k$ nodes, if $-99 < k < 0$ cuts are generated every $-k$ nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if $k=-99$ generate cuts only at the root node, if $k=0$ or $100$ do not generate cuts.");
    roptions->setOptionExtraInfo("2mir_cuts",119);

    roptions->AddLowerBoundedIntegerOption("flow_cover_cuts",
        "Frequency (in terms of nodes) for generating flow cover cuts in branch-and-cut",
        -100,-5,
        "If $k > 0$, cuts are generated every $k$ nodes, if $-99 < k < 0$ cuts are generated every $-k$ nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if $k=-99$ generate cuts only at the root node, if $k=0$ or $100$ do not generate cuts.");
    roptions->setOptionExtraInfo("flow_cover_cuts",119);
    roptions->AddLowerBoundedIntegerOption("lift_and_project_cuts",
        "Frequency (in terms of nodes) for generating lift-and-project cuts in branch-and-cut",
        -100,0,
        "If $k > 0$, cuts are generated every $k$ nodes, if $-99 < k < 0$ cuts are generated every $-k$ nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if $k=-99$ generate cuts only at the root node, if $k=0$ or $100$ do not generate cuts.");
    roptions->setOptionExtraInfo("lift_and_project_cuts", 119);
    roptions->AddLowerBoundedIntegerOption("reduce_and_split_cuts",
        "Frequency (in terms of nodes) for generating reduce-and-split cuts in branch-and-cut",
        -100,0,
        "If $k > 0$, cuts are generated every $k$ nodes, if $-99 < k < 0$ cuts are generated every $-k$ nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if $k=-99$ generate cuts only at the root node, if $k=0$ or $100$ do not generate cuts.");
    roptions->setOptionExtraInfo("reduce_and_split_cuts", 119);


    roptions->AddLowerBoundedIntegerOption("clique_cuts",
        "Frequency (in terms of nodes) for generating clique cuts in branch-and-cut",
        -100,-5,
        "If $k > 0$, cuts are generated every $k$ nodes, if $-99 < k < 0$ cuts are generated every $-k$ nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if $k=-99$ generate cuts only at the root node, if $k=0$ or $100$ do not generate cuts.");
    roptions->setOptionExtraInfo("clique_cuts", 119);

  }


  /** Add milp cut generators according to options.*/
  void
  BonminSetup::addMilpCutGenerators()
  {

    int freq;

    options_->GetIntegerValue("Gomory_cuts", freq,prefix_.c_str());

    if (freq) {
      CuttingMethod cg;
      cg.frequency = freq;
      CglGomory * gom = new CglGomory;
      cg.cgl = gom;
      gom->setLimitAtRoot(5000);
      gom->setLimit(500);
      gom->setLargestFactorMultiplier(1e-08);
      cg.id = "Mixed Integer Gomory";
      cutGenerators_.push_back(cg);
    }

#if 0
    options_->GetIntegerValue("probing_cuts",freq,prefix_.c_str());
    if (freq) {
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
#endif

    options_->GetIntegerValue("mir_cuts",freq,prefix_.c_str());

    if (freq) {
      CuttingMethod cg;
      cg.frequency = freq;
      CglMixedIntegerRounding2 * mir = new CglMixedIntegerRounding2;
      //CglMixedIntegerRounding2 * mir = new CglMixedIntegerRounding2(1, true, 1);
      cg.cgl = mir;
      cg.id = "Mixed Integer Rounding";
      cutGenerators_.push_back(cg);
    }

    options_->GetIntegerValue("2mir_cuts",freq,prefix_.c_str());

    if (freq) {
      CuttingMethod cg;
      cg.frequency = freq;
      CglTwomir * mir2 = new CglTwomir;
      cg.cgl = mir2;
      cg.id = "2-MIR";
      cutGenerators_.push_back(cg);
    }

    options_->GetIntegerValue("cover_cuts",freq,prefix_.c_str());

    if (freq) {
      CuttingMethod cg;
      cg.frequency = freq;
      CglKnapsackCover * cover = new CglKnapsackCover;
      cg.cgl = cover;
      cg.id = "Cover";
      cutGenerators_.push_back(cg);
    }

    options_->GetIntegerValue("clique_cuts",freq,prefix_.c_str());

    if (freq) {
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

    options_->GetIntegerValue("flow_cover_cuts",freq,prefix_.c_str());

    if (freq) {
      CuttingMethod cg;
      cg.frequency = freq;
      CglFlowCover * flow = new CglFlowCover;
      cg.cgl = flow;
      cg.id = "Flow Covers";
      cutGenerators_.push_back(cg);
    }

    options_->GetIntegerValue("lift_and_project_cuts",freq,prefix_.c_str());

    if (freq) {
      CuttingMethod cg;
      cg.frequency = freq;
      CglLandP * landp = new CglLandP;
      cg.cgl = landp;
      cg.id = "Lift-and-Project";
      cutGenerators_.push_back(cg);
    }

    options_->GetIntegerValue("reduce_and_split_cuts",freq,prefix_.c_str());

    if (freq) {
      CuttingMethod cg;
      cg.frequency = freq;
      CglRedSplit * rands = new CglRedSplit;
      cg.cgl = rands;
      cg.id = "Reduce-and-Split";
      cutGenerators_.push_back(cg);
    }
  }


  void
  BonminSetup::initializeBBB()
  {
    continuousSolver_ = nonlinearSolver_;
    nonlinearSolver_->ignoreFailures();
    OsiBabSolver extraStuff(2);
    continuousSolver_->setAuxiliaryInfo(&extraStuff);

    intParam_[BabSetupBase::SpecialOption] = 16;
    if (!options_->GetIntegerValue("number_before_trust",intParam_[BabSetupBase::MinReliability],prefix_.c_str())) {
      intParam_[BabSetupBase::MinReliability] = 1;
      std::string o_name = prefix_ + "number_before_trust";
      options_->SetIntegerValue(o_name.c_str(),intParam_[BabSetupBase::MinReliability],true,true);
    }
    if (!options_->GetIntegerValue("number_strong_branch",intParam_[BabSetupBase::NumberStrong],prefix_.c_str())) {
      intParam_[BabSetupBase::NumberStrong] = 1000;
      std::string o_name = prefix_ + "number_strong_branch";
      options_->SetIntegerValue(o_name.c_str(),intParam_[BabSetupBase::NumberStrong],true,true);
    }
    int varSelection;
    bool val = options_->GetEnumValue("variable_selection",varSelection,prefix_.c_str());
    if (!val){// || varSelection == STRONG_BRANCHING || varSelection == RELIABILITY_BRANCHING ) {
      std::string o_name = prefix_ + "variable_selection";
      options_->SetStringValue(o_name.c_str(), "nlp-strong-branching",true,true);
      varSelection = NLP_STRONG_BRANCHING;
    }

    switch (varSelection) {
#ifdef BONMIN_CURVATURE_BRANCHING
    case CURVATURE_ESTIMATOR:
#endif
    case QP_STRONG_BRANCHING:
    case LP_STRONG_BRANCHING:
    case NLP_STRONG_BRANCHING: {
        continuousSolver_->findIntegersAndSOS(false);
        setPriorities();
        addSos();
        Ipopt::SmartPtr<StrongBranchingSolver> strong_solver = NULL;
        BonChooseVariable * chooseVariable = new BonChooseVariable(*this, nonlinearSolver_);
        chooseVariable->passInMessageHandler(nonlinearSolver_->messageHandler());
        switch (varSelection) {
#ifdef BONMIN_CURVATURE_BRANCHING
        case CURVATURE_ESTIMATOR:
          strong_solver = new CurvBranchingSolver(nonlinearSolver_);
          chooseVariable->setTrustStrongForSolution(false);
          chooseVariable->setTrustStrongForBound(false);
          //chooseVariable->setOnlyPseudoWhenTrusted(true);
          chooseVariable->setOnlyPseudoWhenTrusted(false);
          break;
#endif
        case QP_STRONG_BRANCHING:
          chooseVariable->setTrustStrongForSolution(false);
          strong_solver = new QpBranchingSolver(nonlinearSolver_);
          // The bound returned from the QP can be wrong, since the
          // objective is not guaranteed to be an underestimator:
          chooseVariable->setTrustStrongForBound(false);
          //chooseVariable->setOnlyPseudoWhenTrusted(true);
          chooseVariable->setOnlyPseudoWhenTrusted(false);
          break;
        case LP_STRONG_BRANCHING:
          chooseVariable->setTrustStrongForSolution(false);
          strong_solver = new LpBranchingSolver(this);
          //chooseVariable->setOnlyPseudoWhenTrusted(true);
          chooseVariable->setOnlyPseudoWhenTrusted(false);
          break;
         case NLP_STRONG_BRANCHING:
          chooseVariable->setTrustStrongForSolution(false);
          chooseVariable->setTrustStrongForBound(true);
          chooseVariable->setOnlyPseudoWhenTrusted(false);
          break;
        }
        nonlinearSolver_->SetStrongBrachingSolver(strong_solver);
        branchingMethod_ = chooseVariable;
      }
      break;
    case OSI_SIMPLE:
      continuousSolver_->findIntegersAndSOS(false);
      setPriorities();
      addSos();
      branchingMethod_ = new OsiChooseVariable(nonlinearSolver_);

      break;
    case OSI_STRONG:
      continuousSolver_->findIntegersAndSOS(false);
      setPriorities();
      addSos();
      branchingMethod_ = new OsiChooseStrong(nonlinearSolver_);
      break;
    case RANDOM:
      continuousSolver_->findIntegersAndSOS(false);
      setPriorities();
      addSos();
      branchingMethod_ = new BonRandomChoice(nonlinearSolver_);
      break;
      //default:
      //abort();
    }
    if (branchingMethod_ != NULL) {
      branchingMethod_->setNumberStrong(intParam_[NumberStrong]);
    }


    Ipopt::Index doHeuristicDiveFractional = false;
    options()->GetEnumValue("heuristic_dive_fractional",doHeuristicDiveFractional,prefix_.c_str());
    if(doHeuristicDiveFractional){
      HeuristicDiveFractional* dive_fractional = new HeuristicDiveFractional(this);
      HeuristicMethod h;
      h.heuristic = dive_fractional;
      h.id = "DiveFractional";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicDiveVectorLength = false;
    options()->GetEnumValue("heuristic_dive_vectorLength",doHeuristicDiveVectorLength,prefix_.c_str());
    if(doHeuristicDiveVectorLength){
      HeuristicDiveVectorLength* dive_vectorLength = new HeuristicDiveVectorLength(this);
      HeuristicMethod h;
      h.heuristic = dive_vectorLength;
      h.id = "DiveVectorLength";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicDiveMIPFractional = false;
    if(!options()->GetEnumValue("heuristic_dive_MIP_fractional",doHeuristicDiveMIPFractional,prefix_.c_str())){
      doHeuristicDiveMIPFractional = true;
      std::string o_name = prefix_ + "heuristic_dive_MIP_fractional";
      options_->SetStringValue(o_name.c_str(), "yes",true,true);
    }
    if(doHeuristicDiveMIPFractional){
      HeuristicDiveMIPFractional* dive_MIP_fractional = new HeuristicDiveMIPFractional(this);
      HeuristicMethod h;
      h.heuristic = dive_MIP_fractional;
      h.id = "DiveMIPFractional";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicDiveMIPVectorLength = false;
    options()->GetEnumValue("heuristic_dive_MIP_vectorLength",doHeuristicDiveMIPVectorLength,prefix_.c_str());
    if(doHeuristicDiveMIPVectorLength){
      HeuristicDiveMIPVectorLength* dive_MIP_vectorLength = new HeuristicDiveMIPVectorLength(this);
      HeuristicMethod h;
      h.heuristic = dive_MIP_vectorLength;
      h.id = "DiveMIPVectorLength";
      heuristics_.push_back(h);
    }
    Ipopt::Index doHeuristicFPump = false;
    if(!nonlinearSolver_->model()->hasGeneralInteger() && !options()->GetEnumValue("heuristic_feasibility_pump",doHeuristicFPump,prefix_.c_str())){
      doHeuristicFPump = true;
      std::string o_name = prefix_ + "heuristic_feasibility_pump";
      options_->SetStringValue(o_name.c_str(), "yes",true,true);
    }
    if(doHeuristicFPump){
      HeuristicFPump* feasibility_pump = new HeuristicFPump(this);
      HeuristicMethod h;
      h.heuristic = feasibility_pump;
      h.id = "FPump";
      heuristics_.push_back(h);
    }

    Ipopt::Index doFixAndSolve = false;
    options()->GetEnumValue("fix_and_solve_heuristic",doFixAndSolve,prefix_.c_str());
    if(doFixAndSolve){
      FixAndSolveHeuristic* fix_and_solve = new FixAndSolveHeuristic(this);
      HeuristicMethod h;
      h.heuristic = fix_and_solve;
      h.id = "Fix and Solve";
      heuristics_.push_back(h);
    }

    Ipopt::Index doDummyPump = false;
    options()->GetEnumValue("dummy_pump_heuristic",doDummyPump,prefix_.c_str());
    if(doDummyPump){
      DummyPump* fix_and_solve = new DummyPump(this);
      HeuristicMethod h;
      h.heuristic = fix_and_solve;
      h.id = "Dummy pump";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicRINS = false;
    options()->GetEnumValue("heuristic_RINS",doHeuristicRINS,prefix_.c_str());
    if(doHeuristicRINS){
      HeuristicRINS* rins = new HeuristicRINS(this);
      HeuristicMethod h;
      h.heuristic = rins;
      h.id = "RINS";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicLocalBranching = false;
    options()->GetEnumValue("heuristic_local_branching",doHeuristicLocalBranching,prefix_.c_str());
    if(doHeuristicLocalBranching){
      HeuristicLocalBranching* local_branching = new HeuristicLocalBranching(this);
      HeuristicMethod h;
      h.heuristic = local_branching;
      h.id = "LocalBranching";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicPumpForMinlp = false;
    options()->GetEnumValue("pump_for_minlp",doHeuristicPumpForMinlp,prefix_.c_str());
    if(doHeuristicPumpForMinlp){
      PumpForMinlp * pump = new PumpForMinlp(this);
      HeuristicMethod h;
      h.heuristic = pump;
      h.id = "Pump for MINLP";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicMilpRounding = false;
    options()->GetEnumValue("MILP_rounding_heuristic",doHeuristicMilpRounding,prefix_.c_str());
    if(doHeuristicMilpRounding){
      MilpRounding * round = new MilpRounding(this);
      HeuristicMethod h;
      h.heuristic = round;
      h.id = "MILP Rounding";
      heuristics_.push_back(h);
    }
  }


  void
  BonminSetup::initializeBHyb(bool createContinuousSolver /*= false*/)
  {
    double setup_time = -CoinCpuTime();
    if (createContinuousSolver) {
      /* Create linear solver */
      continuousSolver_ = new OsiClpSolverInterface;
      int lpLogLevel;
      options_->GetIntegerValue("lp_log_level",lpLogLevel,prefix_.c_str());
      if(messageHandler_)
        continuousSolver_->passInMessageHandler(messageHandler_);
      continuousSolver_->messageHandler()->setLogLevel(lpLogLevel);
      nonlinearSolver_->forceSolverOutput(intParam_[RootLogLevel]); 

      if(IsValid(linearizer_))//Use user provided linearizer
        nonlinearSolver_->set_linearizer(linearizer_);

      nonlinearSolver_->extractLinearRelaxation(*continuousSolver_);
      nonlinearSolver_->setSolverOutputToDefault(); 

      // say bound dubious, does cuts at solution
      OsiBabSolver * extraStuff = new OsiBabSolver(3);
      continuousSolver_->setAuxiliaryInfo(extraStuff);
      delete extraStuff;
    }
    Algorithm algo = getAlgorithm();
    std::string prefix = (prefix_ == "bonmin.") ? "" : prefix_;
    if (algo == B_Hyb) {
      std::string o_name = prefix_ + "oa_decomposition";
      options_->SetStringValue(o_name.c_str(),"no", true, true);
      o_name = prefix_ + "pump_for_minlp";
      options_->SetStringValue(o_name.c_str(),"yes", true, true);
      o_name = prefix + "pump_for_minlp.time_limit";
      options_->SetNumericValue(o_name.c_str(),10, true, true);
      o_name = prefix + "pump_for_minlp.solution_limit";
      options_->SetIntegerValue(o_name.c_str(),3, true, true);
    }
    else if (algo == B_OA) {
      std::string o_name = prefix_ + "oa_decomposition";
      options_->SetStringValue(o_name.c_str(),"yes", true, true);
      o_name = prefix + "oa_decomposition.time_limit";
      options_->SetNumericValue(o_name.c_str(),DBL_MAX, true, true);
      o_name = prefix_ + "pump_for_minlp";
      options_->SetStringValue(o_name.c_str(),"no", true, true);
      o_name = prefix + "nlp_solve_frequency";
      options_->SetIntegerValue(o_name.c_str(), 0, true, true);
      o_name = prefix + "bb_log_level";
      options_->SetIntegerValue(o_name.c_str(), 0, true, true);
    }
    else if (algo == B_IFP) {
      std::string o_name = prefix_ + "oa_decomposition";
      options_->SetStringValue(o_name.c_str(),"no", true, true);
      o_name = prefix_ + "pump_for_minlp";
      options_->SetStringValue(o_name.c_str(),"yes", true, true);
      o_name = prefix + "pump_for_minlp.time_limit";
      options_->SetNumericValue(o_name.c_str(),DBL_MAX, true, true);
      o_name = prefix_ + "nlp_solve_frequency";
      options_->SetIntegerValue(o_name.c_str(), 0, true, true);
      o_name = prefix_ + "fp_pass_infeasible";
      options_->SetStringValue(o_name.c_str(), "yes", true, true);
      //o_name = prefix_ + "cutoff_decr";
      //options_->SetNumericValue(o_name.c_str(), 1e-02, true, true);
      intParam_[BabLogLevel] = 0;
    }
    else if (algo==B_QG) {
      std::string o_name = prefix_ + "oa_decomposition";
      options_->SetStringValue(o_name.c_str(),"no", true, true);
      o_name = prefix_ + "pump_for_minlp";
      options_->SetStringValue(o_name.c_str(),"no", true, true);
      o_name = prefix_ + "nlp_solve_frequency";
      options_->SetIntegerValue(o_name.c_str(), 0, true, true);
    }
    else if (algo==B_Ecp) {
      std::string o_name = prefix_ + "oa_decomposition";
      options_->SetStringValue(o_name.c_str(),"no", true, true);
      o_name = prefix_ + "pump_for_minlp";
      options_->SetStringValue(o_name.c_str(),"no", true, true);
      o_name = prefix_ + "nlp_solve_frequency";
      options_->SetIntegerValue(o_name.c_str(), 0, true, true);
      o_name = prefix_ + "filmint_ecp_cuts";
      options_->SetIntegerValue(o_name.c_str(), 1, true, true);
    }
//#define GREAT_STUFF_FOR_ANDREAS
#ifdef GREAT_STUFF_FOR_ANDREAS
    printf("ToDo: Clean me up in Bab::branchAndBound\n");
    OsiCuts cuts;
    nonlinearSolver_->getOuterApproximation(cuts, true, NULL, true);
    continuousSolver_->applyCuts(cuts);
#endif


    int varSelection;
    options_->GetEnumValue("variable_selection",varSelection,prefix_.c_str());
    if (varSelection > RELIABILITY_BRANCHING) {
      switch (varSelection){
        case OSI_SIMPLE:
          continuousSolver_->findIntegersAndSOS(false);
          setPriorities();
          addSos();
          branchingMethod_ = new OsiChooseVariable(nonlinearSolver_);
    
          break;
        case OSI_STRONG:
          {
          continuousSolver_->findIntegersAndSOS(false);
          setPriorities();
          addSos();
          OsiChooseStrong * chooser = new OsiChooseStrong(nonlinearSolver_);
          branchingMethod_ = chooser;
          chooser->setNumberStrong(intParam_[NumberStrong]);
          chooser->setTrustStrongForSolution(false);
          chooser->setNumberBeforeTrusted(intParam_[MinReliability]);
          }
          break;
        default:
          std::cout<<"Variable selection stragey not available with oa branch-and-cut."<<std::endl;
          break;
     }
    }
    /* Populate cut generation and heuristic procedures.*/
    int ival;
    options_->GetIntegerValue("nlp_solve_frequency",ival,prefix_.c_str());
    if (ival != 0) {
      CuttingMethod cg;
      cg.frequency = ival;
      OaNlpOptim * nlpsol = new OaNlpOptim(*this);
      nlpsol->passInMessageHandler(messageHandler_);
      cg.cgl = nlpsol;
      cg.id="NLP solution based oa cuts";
      cutGenerators_.push_back(cg);
    }

    options_->GetIntegerValue("filmint_ecp_cuts",ival, prefix_.c_str());
    if (ival != 0) {
      CuttingMethod cg;
      cg.frequency = ival;
      EcpCuts * ecp = new EcpCuts(*this);
      ecp->passInMessageHandler(messageHandler_);
      cg.cgl = ecp;
      cg.id = "Ecp cuts";
      cutGenerators_.push_back(cg);
    }

    if (algo == B_Hyb || algo == B_Ecp)
      addMilpCutGenerators();

    int doFp;
    options_->GetEnumValue("pump_for_minlp",doFp,prefix_.c_str());
    if (doFp) {
      CuttingMethod cg;
      cg.frequency = -99;
      MinlpFeasPump * oa = new MinlpFeasPump(*this);
      oa->passInMessageHandler(messageHandler_);
      cg.cgl = oa;
      cg.id = "Feasibility Pump for MINLP.";
      cutGenerators_.push_back(cg);

    }
    int doOa;
    options_->GetEnumValue("oa_decomposition",doOa,prefix_.c_str());
    if (doOa) {
      CuttingMethod cg;
      cg.frequency = -99;
      OACutGenerator2 * oa = new OACutGenerator2(*this);
      oa->passInMessageHandler(messageHandler_);
      cg.cgl = oa;
      cg.id = "Outer Approximation decomposition.";
      cutGenerators_.push_back(cg);

    }

    {
      CuttingMethod cg;
      cg.frequency = 1;
      OaFeasibilityChecker * oa = new OaFeasibilityChecker(*this);
      oa->passInMessageHandler(messageHandler_);
      oa->setReassignLpSolver(false);
      cg.cgl = oa;
      cg.id = "Outer Approximation feasibility check.";
      cg.atSolution = false;
      cg.normal = true;
      cg.always = true;
      cutGenerators_.push_back(cg);
    }

    {
      CuttingMethod cg;
      cg.frequency = 1;
      OaFeasibilityChecker * oa = new OaFeasibilityChecker(*this);
      oa->passInMessageHandler(messageHandler_);
      oa->setReassignLpSolver(true);
      cg.cgl = oa;
      cg.id = "Outer Approximation strong branching solution check.";
      cg.atSolution = true;
      cg.normal = false;
      cutGenerators_.push_back(cg);
    }

    DummyHeuristic * oaHeu = new DummyHeuristic;
    oaHeu->setNlp(nonlinearSolver_);
    HeuristicMethod h;
    h.heuristic = oaHeu;
    h.id = "nonlinear program";
    heuristics_.push_back(h);

    Ipopt::Index doHeuristicRINS = false;
    options()->GetEnumValue("heuristic_RINS",doHeuristicRINS,prefix_.c_str());
    if(doHeuristicRINS){
      HeuristicRINS* rins = new HeuristicRINS(this);
      HeuristicMethod h;
      h.heuristic = rins;
      h.id = "RINS";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicLocalBranching = false;
    options()->GetEnumValue("heuristic_local_branching",doHeuristicLocalBranching,prefix_.c_str());
    if(doHeuristicLocalBranching){
      HeuristicLocalBranching* local_branching = new HeuristicLocalBranching(this);
      HeuristicMethod h;
      h.heuristic = local_branching;
      h.id = "LocalBranching";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicFPump = false;
    options()->GetEnumValue("heuristic_feasibility_pump",doHeuristicFPump,prefix_.c_str());
    if(doHeuristicFPump){
      HeuristicFPump* feasibility_pump = new HeuristicFPump(this);
      HeuristicMethod h;
      h.heuristic = feasibility_pump;
      h.id = "FPump";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicDiveFractional = false;
    options()->GetEnumValue("heuristic_dive_fractional",doHeuristicDiveFractional,prefix_.c_str());
    if(doHeuristicDiveFractional){
      HeuristicDiveFractional* dive_fractional = new HeuristicDiveFractional(this);
      HeuristicMethod h;
      h.heuristic = dive_fractional;
      h.id = "DiveFractional";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicDiveVectorLength = false;
    options()->GetEnumValue("heuristic_dive_vectorLength",doHeuristicDiveVectorLength,prefix_.c_str());
    if(doHeuristicDiveVectorLength){
      HeuristicDiveVectorLength* dive_vectorLength = new HeuristicDiveVectorLength(this);
      HeuristicMethod h;
      h.heuristic = dive_vectorLength;
      h.id = "DiveVectorLength";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicDiveMIPFractional = false;
    options()->GetEnumValue("heuristic_dive_MIP_fractional",doHeuristicDiveMIPFractional,prefix_.c_str());
    if(doHeuristicDiveMIPFractional){
      HeuristicDiveMIPFractional* dive_MIP_fractional = new HeuristicDiveMIPFractional(this);
      HeuristicMethod h;
      h.heuristic = dive_MIP_fractional;
      h.id = "DiveMIPFractional";
      heuristics_.push_back(h);
    }

    Ipopt::Index doHeuristicDiveMIPVectorLength = false;
    options()->GetEnumValue("heuristic_dive_MIP_vectorLength",doHeuristicDiveMIPVectorLength,prefix_.c_str());
    if(doHeuristicDiveMIPVectorLength){
      HeuristicDiveMIPVectorLength* dive_MIP_vectorLength = new HeuristicDiveMIPVectorLength(this);
      HeuristicMethod h;
      h.heuristic = dive_MIP_vectorLength;
      h.id = "DiveMIPVectorLength";
      heuristics_.push_back(h);
    }

#if 0
    if(true){
      InnerApproximation * inner = new InnerApproximation(this);
      HeuristicMethod h;
      h.heuristic = inner;
      h.id = "InnerApproximation";
      heuristics_.push_back(h);
    }
#endif
    setup_time += CoinCpuTime();
    doubleParam_[MaxTime] -= setup_time;
  }


  Algorithm BonminSetup::getAlgorithm()
  {
    if (algo_ != Dummy)
      return algo_;
    if (IsValid(options_)) {
      int ival;
      options_->GetEnumValue("algorithm", ival,prefix_.c_str());
      return Algorithm(ival);
    }
    else return Algorithm(3);
  }

}/* end namespace Bonmin*/


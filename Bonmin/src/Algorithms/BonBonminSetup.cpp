// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/13/2007

#include "BonminConfig.h"
#include "OsiClpSolverInterface.hpp"

#include "BonBonminSetup.hpp"
#include "BonCurvBranchingSolver.hpp"
#include "BonChooseVariable.hpp"
#include "BonRandomChoice.hpp"
#include "BonDiver.hpp"
#include "BonQpBranchingSolver.hpp"
#include "BonLpBranchingSolver.hpp"

//OA machinery
#include "BonDummyHeuristic.hpp"
#include "BonOACutGenerator2.hpp"
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

namespace Bonmin
{
  BonminSetup::BonminSetup():BabSetupBase(),algo_(Dummy)
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
      options_->GetIntegerValue("lp_log_level",lpLogLevel,"bonmin.");
      lpMessageHandler_ = nonlinearSolver_->messageHandler()->clone();
      continuousSolver_->passInMessageHandler(lpMessageHandler_);
      continuousSolver_->messageHandler()->setLogLevel(lpLogLevel);
      nonlinearSolver_->extractLinearRelaxation(*continuousSolver_);
      // say bound dubious, does cuts at solution
      OsiBabSolver * extraStuff = new OsiBabSolver(3);
      continuousSolver_->setAuxiliaryInfo(extraStuff);
      delete extraStuff;
    }
  }
  void BonminSetup::registerAllOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    BabSetupBase::registerAllOptions(roptions);

    /* Outer Approximation options.*/
    OACutGenerator2::registerOptions(roptions);
    EcpCuts::registerOptions(roptions);
    OaNlpOptim::registerOptions(roptions);


    BonCbcFullNodeInfo::registerOptions(roptions);


    registerMilpCutGenerators(roptions);


    roptions->SetRegisteringCategory("Algorithm choice", RegisteredOptions::BonminCategory);
    roptions->AddStringOption5("algorithm",
        "Choice of the algorithm.",
        "B-BB",
        "B-BB","simple branch-and-bound algorithm,",
        "B-OA","OA Decomposition algorithm,",
        "B-QG","Quesada and Grossmann branch-and-cut algorithm,",
        "B-Hyb","hybrid outer approximation based branch-and-cut,",
        "B-Ecp","ecp cuts based branch-and-cut a la FilMINT.",
        "This will preset some of the options of bonmin depending on the algorithm choice."
                              );
    roptions->setOptionExtraInfo("algorithm",31);
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
    roptions->SetRegisteringCategory("MILP cutting planes in hybrid", RegisteredOptions::BonminCategory);

    roptions->AddLowerBoundedIntegerOption("Gomory_cuts",
        "Frequency k (in terms of nodes) for generating Gomory cuts in branch-and-cut.",
        -100,-5,
        "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
    roptions->setOptionExtraInfo("Gomory_cuts",5);
#if 0
    roptions->AddBoundedIntegerOption("probing_cuts",
        "Frequency (in terms of nodes) for generating probing cuts in branch-and-cut",
        0,0,0,
        "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
    roptions->setOptionExtraInfo("probing_cuts",5);
#endif 
    roptions->AddLowerBoundedIntegerOption("cover_cuts",
        "Frequency (in terms of nodes) for generating cover cuts in branch-and-cut",
        -100,-5,
        "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
    roptions->setOptionExtraInfo("cover_cuts",5);

    roptions->AddLowerBoundedIntegerOption("mir_cuts",
        "Frequency (in terms of nodes) for generating MIR cuts in branch-and-cut",
        -100,-5,
        "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
    roptions->setOptionExtraInfo("mir_cuts",5);
    roptions->AddLowerBoundedIntegerOption("2mir_cuts",
        "Frequency (in terms of nodes) for generating 2-MIR cuts in branch-and-cut",
        -100,0,
        "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
    roptions->setOptionExtraInfo("2mir_cuts",5);
    roptions->AddLowerBoundedIntegerOption("flow_covers_cuts",
        "Frequency (in terms of nodes) for generating flow cover cuts in branch-and-cut",
        -100,-5,
        "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
    roptions->setOptionExtraInfo("flow_covers_cuts",5);
    roptions->AddLowerBoundedIntegerOption("lift_and_project_cuts",
        "Frequency (in terms of nodes) for generating lift-and-project cuts in branch-and-cut",
        -100,0,
        "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
    roptions->setOptionExtraInfo("lift_and_project_cuts",5);
    roptions->AddLowerBoundedIntegerOption("reduce_and_split_cuts",
        "Frequency (in terms of nodes) for generating reduce-and-split cuts in branch-and-cut",
        -100,0,
        "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
    roptions->setOptionExtraInfo("reduce_and_split_cuts",5);
    roptions->AddLowerBoundedIntegerOption("clique_cuts",
        "Frequency (in terms of nodes) for generating clique cuts in branch-and-cut",
        -100,-5,
        "If k > 0, cuts are generated every k nodes, if -99 < k < 0 cuts are generated every -k nodes but "
        "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
        "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");
    roptions->setOptionExtraInfo("clique_cuts",5);
  }
  /** Add milp cut generators according to options.*/
  void
  BonminSetup::addMilpCutGenerators()
  {
    int freq;
    options_->GetIntegerValue("Gomory_cuts", freq,"bonmin.");
    if (freq) {
      CuttingMethod cg;
      cg.frequency = freq;
      CglGomory * gom = new CglGomory;
      cg.cgl = gom;
      gom->setLimitAtRoot(512);
      gom->setLimit(50);
      cg.id = "Mixed Integer Gomory";
      cutGenerators_.push_back(cg);
    }
#if 0
    options_->GetIntegerValue("probing_cuts",freq,"bonmin.");
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
    options_->GetIntegerValue("mir_cuts",freq,"bonmin.");
    if (freq) {
      CuttingMethod cg;
      cg.frequency = freq;
      CglMixedIntegerRounding2 * mir = new CglMixedIntegerRounding2;
      cg.cgl = mir;
      cg.id = "Mixed Integer Rounding";
      cutGenerators_.push_back(cg);


    }
    options_->GetIntegerValue("2mir_cuts",freq,"bonmin.");
    if (freq) {
      CuttingMethod cg;
      cg.frequency = freq;
      CglTwomir * mir2 = new CglTwomir;
      cg.cgl = mir2;
      cg.id = "2-MIR";
      cutGenerators_.push_back(cg);
    }
    options_->GetIntegerValue("cover_cuts",freq,"bonmin.");
    if (freq) {
      CuttingMethod cg;
      cg.frequency = freq;
      CglKnapsackCover * cover = new CglKnapsackCover;
      cg.cgl = cover;
      cg.id = "Cover";
      cutGenerators_.push_back(cg);
    }

    options_->GetIntegerValue("clique_cuts",freq,"bonmin.");
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
    options_->GetIntegerValue("flow_covers_cuts",freq,"bonmin.");
    if (freq) {
      CuttingMethod cg;
      cg.frequency = freq;
      CglFlowCover * flow = new CglFlowCover;
      cg.cgl = flow;
      cg.id = "Flow Covers";
      cutGenerators_.push_back(cg);
    }
    options_->GetIntegerValue("lift_and_project_cuts",freq,"bonmin.");
    if (freq) {
      CuttingMethod cg;
      cg.frequency = freq;
      CglLandP * landp = new CglLandP;
      cg.cgl = landp;
      cg.id = "Lift-and-Project";
      cutGenerators_.push_back(cg);
    }
    options_->GetIntegerValue("reduce_and_split_cuts",freq,"bonmin.");
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
#if 1
    if (!options_->GetIntegerValue("number_before_trust",intParam_[BabSetupBase::MinReliability],"bonmin.")) {
      intParam_[BabSetupBase::MinReliability] = 1;
      options_->SetIntegerValue("bonmin.number_before_trust",intParam_[BabSetupBase::MinReliability], true, true);
    }
    if (!options_->GetIntegerValue("number_strong_branch",intParam_[BabSetupBase::NumberStrong],"bonmin.")) {
      intParam_[BabSetupBase::NumberStrong] = 1000;
      options_->SetIntegerValue("bonmin.number_strong_branch",intParam_[BabSetupBase::NumberStrong], true, true);
    }
    int varSelection;
    bool val = options_->GetEnumValue("variable_selection",varSelection,"bonmin.");
    // Set branching strategy
    if (varSelection == MOST_FRACTIONAL) {
      intParam_[NumberStrong] = 0;
      intParam_[MinReliability] = 0;
      options_->SetIntegerValue("bonmin.number_strong_branch",intParam_[BabSetupBase::NumberStrong],true, true);
    }
    if (!val || varSelection == STRONG_BRANCHING || varSelection == RELIABILITY_BRANCHING ) {
      options_->SetStringValue("bonmin.variable_selection", "nlp-strong-branching", true, true);
      varSelection = NLP_STRONG_BRANCHING;
    }
#endif

    switch (varSelection) {
    case CURVATURE_ESTIMATOR:
    case QP_STRONG_BRANCHING:
    case LP_STRONG_BRANCHING:
    case NLP_STRONG_BRANCHING: {
        continuousSolver_->findIntegersAndSOS(false);
        setPriorities();
        addSos();
        SmartPtr<StrongBranchingSolver> strong_solver = NULL;
        BonChooseVariable * chooseVariable = new BonChooseVariable(*this, nonlinearSolver_);
        chooseVariable->passInMessageHandler(nonlinearSolver_->messageHandler());
        switch (varSelection) {
        case CURVATURE_ESTIMATOR:
          strong_solver = new CurvBranchingSolver(nonlinearSolver_);
          chooseVariable->setTrustStrongForSolution(false);
          chooseVariable->setTrustStrongForBound(false);
          //chooseVariable->setOnlyPseudoWhenTrusted(true);
          chooseVariable->setOnlyPseudoWhenTrusted(false);
          break;
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
          strong_solver = new LpBranchingSolver(nonlinearSolver_);
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
  }

  void
  BonminSetup::initializeBHyb(bool createContinuousSolver /*= false*/)
  {
    if (createContinuousSolver) {
      /* Create linear solver */
      continuousSolver_ = new OsiClpSolverInterface;
      int lpLogLevel;
      options_->GetIntegerValue("lp_log_level",lpLogLevel,"bonmin.");
      lpMessageHandler_ = nonlinearSolver_->messageHandler()->clone();
      continuousSolver_->passInMessageHandler(lpMessageHandler_);
      continuousSolver_->messageHandler()->setLogLevel(lpLogLevel);
      nonlinearSolver_->extractLinearRelaxation(*continuousSolver_);
      // say bound dubious, does cuts at solution
      OsiBabSolver * extraStuff = new OsiBabSolver(3);
      continuousSolver_->setAuxiliaryInfo(extraStuff);
      delete extraStuff;
    }
    Algorithm algo = getAlgorithm();
    if (algo == B_OA) {
      options_->SetNumericValue("bonmin.oa_dec_time_limit",COIN_DBL_MAX, true, true);
      options_->SetIntegerValue("bonmin.nlp_solve_frequency", 0, true, true);
      intParam_[BabLogLevel] = 0;
    }
    else if (algo==B_QG) {
      options_->SetNumericValue("bonmin.oa_dec_time_limit",0, true, true);
      options_->SetIntegerValue("bonmin.nlp_solve_frequency", 0, true, true);
    }
    else if (algo==B_Ecp) {
      options_->SetNumericValue("bonmin.oa_dec_time_limit",0, true, true);
      options_->SetIntegerValue("bonmin.nlp_solve_frequency", 0, true, true);
      options_->SetIntegerValue("bonmin.filmint_ecp_cuts", 1, true, true);
      options_->SetIntegerValue("bonmin.number_cut_passes", 1, true, true);
    }
//#define GREAT_STUFF_FOR_ANDREAS
#ifdef GREAT_STUFF_FOR_ANDREAS
    printf("ToDo: Clean me up in Bab::branchAndBound\n");
    OsiCuts cuts;
    nonlinearSolver_->getOuterApproximation(cuts, true, NULL, true);
    continuousSolver_->applyCuts(cuts);
#endif


    int varSelection;
    options_->GetEnumValue("variable_selection",varSelection,"bonmin.");
    if (varSelection > RELIABILITY_BRANCHING) {
      std::cout<<"Variable selection stragey not available with oa branch-and-cut."<<std::endl;
    }
    /* Populate cut generation and heuristic procedures.*/
    int ival;
    options_->GetIntegerValue("nlp_solve_frequency",ival,"bonmin.");
    if (ival != 0) {
      CuttingMethod cg;
      cg.frequency = ival;
      OaNlpOptim * nlpsol = new OaNlpOptim(*this);
      nlpsol->passInMessageHandler(nonlinearSolver_->messageHandler());
      cg.cgl = nlpsol;
      cg.id="NLP solution based oa cuts";
      cutGenerators_.push_back(cg);
    }

    options_->GetIntegerValue("filmint_ecp_cuts",ival, "bonmin.");
    if (ival != 0) {
      CuttingMethod cg;
      cg.frequency = ival;
      EcpCuts * ecp = new EcpCuts(*this);
      ecp->passInMessageHandler(nonlinearSolver_->messageHandler());
      cg.cgl = ecp;
      cg.id = "Ecp cuts";
      cutGenerators_.push_back(cg);
    }

    if (algo == B_Hyb || algo == B_Ecp)
      addMilpCutGenerators();

    double oaTime;
    options_->GetNumericValue("oa_dec_time_limit",oaTime,"bonmin.");
    if (oaTime > 0.) {
      CuttingMethod cg;
      cg.frequency = -99;
      OACutGenerator2 * oa = new OACutGenerator2(*this);
      oa->passInMessageHandler(nonlinearSolver_->messageHandler());
      cg.cgl = oa;
      cg.id = "Outer Approximation decomposition.";
      cutGenerators_.push_back(cg);

    }

    {
      CuttingMethod cg;
      cg.frequency = 1;
      OaFeasibilityChecker * oa = new OaFeasibilityChecker(*this);
      oa->passInMessageHandler(nonlinearSolver_->messageHandler());
      cg.cgl = oa;
      cg.id = "Outer Approximation feasibility check.";
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
  }

  Algorithm BonminSetup::getAlgorithm()
  {
    if (algo_ != Dummy)
      return algo_;
    if (IsValid(options_)) {
      int ival;
      options_->GetEnumValue("algorithm", ival,"bonmin.");
      return Algorithm(ival);
    }
    else return Algorithm(3);
  }

}/* end namespace Bonmin*/


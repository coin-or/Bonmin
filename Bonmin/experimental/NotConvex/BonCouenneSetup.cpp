// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/18/2007

#include "BonCouenneSetup.hpp"
#include "BonNlpHeuristic.hpp"
#include "BonInitHeuristic.hpp"
#include "BonCouenneInterface.hpp"

#include "CouenneObject.hpp"
#include "CouenneVarObject.hpp"
#include "CouenneVTObject.hpp"
#include "CouenneChooseVariable.hpp"
#include "CouenneChooseStrong.hpp"
#include "CouenneSolverInterface.hpp"
#include "CouenneCutGenerator.hpp"
#include "BonCouenneInfo.hpp"
#include "BonCbcNode.hpp"

// MILP cuts
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

// Ampl includes
#include "asl.h"
#include "getstub.h"


namespace Bonmin{
  
  SmartAsl::~SmartAsl(){
    //Code from Ipopt::AmplTNLP to free asl
    if(asl != NULL){
        if (X0) {
          delete [] X0;
          X0 = NULL;
        }
        if (havex0) {
          delete [] havex0;
          havex0 = NULL;
        }
        if (pi0) {
          delete [] pi0;
          pi0 = NULL;
        }
        if (havepi0) {
          delete [] havepi0;
          havepi0 = NULL;
        }
        ASL* asl_to_free = (ASL*)asl;
        ASL_free(&asl_to_free);
        asl = NULL;
    }
    ASL_free(&asl);
  }
  
  CouenneSetup::~CouenneSetup(){
    //if (CouennePtr_)
    //delete CouennePtr_;
  }

  void CouenneSetup::InitializeCouenne (char **& argv) {
    /* Get the basic options. */
    readOptionsFile();
 
    /** Change default value for failure behavior so that code doesn't crash 
	when Ipopt does not solve a sub-problem.*/

    options_->SetStringValue ("nlp_failure_behavior", "fathom", "bonmin.");

    gatherParametersValues(options_);

    continuousSolver_ = new CouenneSolverInterface;
    CouenneInterface * ci = new CouenneInterface;
    nonlinearSolver_ = ci;
    /* Read the model in various places. */
    ci->readAmplNlFile(argv,roptions(),options(),journalist());
    aslfg_ = new SmartAsl;
    aslfg_->asl = readASLfg (argv);

    /** Set the output level for the journalist for all Couenne
     categories.  We probably want to make that a bit more flexible
     later. */
    int i;

    options()->GetIntegerValue("boundtightening_print_level", i, "bonmin.");
    journalist()->GetJournal("console")->
      SetPrintLevel(J_BOUNDTIGHTENING, (EJournalLevel)i);
    options()->GetIntegerValue("branching_print_level", i, "bonmin.");
    journalist()->GetJournal("console")->
      SetPrintLevel(J_BRANCHING, (EJournalLevel)i);
    options()->GetIntegerValue("convexifying_print_level", i, "bonmin.");
    journalist()->GetJournal("console")->
      SetPrintLevel(J_CONVEXIFYING, (EJournalLevel)i);
    options()->GetIntegerValue("problem_print_level", i, "bonmin.");
    journalist()->GetJournal("console")->
      SetPrintLevel(J_PROBLEM, (EJournalLevel)i);

    /* Initialize Couenne cut generator.*/
    //int ivalue, num_points;
    //options()->GetEnumValue("convexification_type", ivalue,"bonmin.");
    //options()->GetIntegerValue("convexification_points",num_points,"bonmin.");
    
    CouenneCutGenerator * couenneCg = 
      new CouenneCutGenerator (ci, this, aslfg_->asl, journalist());

    CouenneProblem * couenneProb = couenneCg -> Problem();

    Bonmin::BabInfo * extraStuff = new Bonmin::CouenneInfo(0);

    // as per instructions by John Forrest, to get changed bounds
    extraStuff -> setExtraCharacteristics (extraStuff -> extraCharacteristics () | 2);

    continuousSolver_ -> setAuxiliaryInfo (extraStuff);
    delete extraStuff;
    
    extraStuff = dynamic_cast <Bonmin::BabInfo *> (continuousSolver_ -> getAuxiliaryInfo ());
    
    /* Setup log level*/
    int lpLogLevel;
    options()->GetIntegerValue("lp_log_level",lpLogLevel,"bonmin.");
    continuousSolver_->messageHandler()->setLogLevel(lpLogLevel);

    //////////////////////////////////////////////////////////////

    ci -> extractLinearRelaxation (*continuousSolver_, *couenneCg);

    // In case there are no discrete variables, we have already a
    // heuristic solution for which create a initialization heuristic
    if (ci -> isProvenOptimal () && 
	ci -> haveNlpSolution ()) {

      /// setup initial heuristic (in principle it should only run once...)
      InitHeuristic* initHeuristic = new InitHeuristic 
	(ci -> getObjValue (), ci -> getColSolution (), *couenneProb);

      heuristics_.push_back(initHeuristic);
    }

    if(extraStuff->infeasibleNode()){
      std::cout<<"Initial linear relaxation constructed by Couenne is infeasible, quit"<<std::endl;
      return;
    }

    continuousSolver_ -> findIntegersAndSOS (false);

    //model -> assignSolver (continuousSolver_, true);
    //continuousSolver_ = model -> solver();

    // add Couenne objects for branching /////////////////////////////////////////////

    std::string s;

    options () -> GetStringValue ("display_stats", s, "couenne.");
    displayStats_ = (s == "yes");

    options () -> GetStringValue ("branching_object", s, "couenne.");

    enum CouenneObject::branch_obj objType = CouenneObject::VAR_OBJ;

    if      (s == "vt_obj")   objType = CouenneObject::VT_OBJ;
    else if (s == "var_obj")  objType = CouenneObject::VAR_OBJ;
    else if (s == "expr_obj") objType = CouenneObject::EXPR_OBJ;
    else {
      printf ("CouenneSetup: Unknown branching object type\n");
      exit (-1);
    }

    int nAuxs = 0, nobj = 0,
      nVars = couenneProb -> nVars ();

    nAuxs = couenneProb -> nVars ();

    OsiObject ** objects = new OsiObject* [nAuxs];

    for (int i = 0; i < nVars; i++) { // for each variable

      exprVar *var = couenneProb -> Var (i);

      // we only want enabled variables
      if (var -> Multiplicity () == 0) 
	continue;

      switch (objType) {

      case CouenneObject::EXPR_OBJ:

	// if this variable is associated with a nonlinear function
	if ((var -> Type  () == AUX) && 
	    (var -> Image () -> Linearity () > LINEAR)) {

	  objects [nobj] = new CouenneObject (couenneProb, var, this, journalist ());
	  objects [nobj++] -> setPriority (1);
	}

	break;

      case CouenneObject::VAR_OBJ:

	// branching objects on variables
	if // comment three lines below for linear variables too
	  (couenneProb -> Dependence () [var -> Index ()] . size () > 0) {  // has indep
	   //|| ((var -> Type () == AUX) &&                                  // or, aux 
	   //    (var -> Image () -> Linearity () > LINEAR))) {              // of nonlinear

	  objects [nobj] = new CouenneVarObject (couenneProb, var, this, journalist ());
	  objects [nobj++] -> setPriority (1);
	}

	break;

      default:
      case CouenneObject::VT_OBJ:

	// branching objects on variables
	if // comment three lines below for linear variables too
	  (couenneProb -> Dependence () [var -> Index ()] . size () > 0) { // has indep
	  //|| ((var -> Type () == AUX) &&                      // or, aux 
	  //(var -> Image () -> Linearity () > LINEAR))) { // of nonlinear

	  objects [nobj] = new CouenneVTObject (couenneProb, var, this, journalist ());
	  objects [nobj++] -> setPriority (1);
	}

	break;
      }
    }

    // Add objects /////////////////////////////////
    continuousSolver_ -> addObjects (nobj, objects);

    for (int i = 0 ; i < nobj ; i++)
      delete objects [i];

    delete [] objects;

    // Setup Convexifier generators ////////////////////////////////////////////////

    int freq;
    options()->GetIntegerValue("convexification_cuts",freq,"couenne.");
    if (freq != 0) {
      CuttingMethod cg;
      cg.frequency = freq;
      cg.cgl = couenneCg;
      cg.id = "Couenne convexifier cuts";
      cutGenerators().push_back(cg);

      // set cut gen pointer
      dynamic_cast <CouenneSolverInterface *> 
	(continuousSolver_) -> setCutGenPtr (couenneCg);
    }

    // add other cut generators -- test for integer variables first
    if (couenneCg -> Problem () -> nIntVars () > 0)
      addMilpCutGenerators ();

    CouennePtr_ = couenneCg;

    /*Setup heuristic to solve nlp problems.*/
    int doNlpHeurisitic = 0;
    options()->GetEnumValue("local_optimization_heuristic", doNlpHeurisitic, "couenne.");
    if(doNlpHeurisitic)
    {
      int numSolve;
      options()->GetIntegerValue("log_num_local_optimization_per_level",numSolve,"couenne.");
      NlpSolveHeuristic * nlpHeuristic = new NlpSolveHeuristic;
      nlpHeuristic->setNlp(*ci,false);
      nlpHeuristic->setCouenneProblem(couenneProb);
      //nlpHeuristic->setMaxNlpInf(1e-4);
      nlpHeuristic->setMaxNlpInf(maxNlpInf_0);
      nlpHeuristic->setNumberSolvePerLevel(numSolve);
      heuristics_.push_back(nlpHeuristic);
    }

    int varSelection;
    if (!options_->GetEnumValue("variable_selection",varSelection,"bonmin.")) {
      // change the default for Couenne
      varSelection = OSI_SIMPLE;
    }

    /*
    if ((varSelection == OSI_STRONG) && 
	(objType == CouenneObject::VT_OBJ)) {
      printf ("Couenne: warning, strong branching is incompatible with Violation Transfer\n"
	      "Resetting variable selection to simple branching");
      varSelection = OSI_SIMPLE;
    }
    */

    switch (varSelection) {

    case OSI_STRONG: { // strong branching
      CouenneChooseStrong * chooseVariable = new CouenneChooseStrong
	(*this, couenneProb, journalist ());
      chooseVariable->setTrustStrongForSolution(false);
      chooseVariable->setTrustStrongForBound(false);
      chooseVariable->setOnlyPseudoWhenTrusted(true);
      branchingMethod_ = chooseVariable;
      break;
    }

    case OSI_SIMPLE: // default choice
      branchingMethod_ = new CouenneChooseVariable 
	(continuousSolver_, couenneProb, journalist ());
      break;

    default:
      std::cerr << "Unknown variable_selection for Couenne\n" << std::endl;
      throw;
      break;
    }

    int ival;
    if (!options_->GetEnumValue("node_comparison",ival,"bonmin.")) {
      // change default for Couenne
      nodeComparisonMethod_ = bestBound;
    }
    else {
      nodeComparisonMethod_ = NodeComparison(ival);
    }

    if(intParam_[NumCutPasses] < 2)
    intParam_[NumCutPasses] = 2;

    // Tell Cbc not to check again if a solution returned from
    // heuristic is indeed feasible
    intParam_[BabSetupBase::SpecialOption] = 16 | 4;
  }
 
  void CouenneSetup::registerOptions(){
    registerAllOptions(roptions());
  }


  void
  CouenneSetup::registerAllOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
    BabSetupBase::registerAllOptions(roptions);
    BonCbcFullNodeInfo::registerOptions(roptions);
    CouenneCutGenerator::registerOptions (roptions);

    roptions -> AddStringOption2 (
      "display_stats",
      "display statistics at the end of the run",
      "no",
      "yes", "",
      "no", "");

    roptions->AddBoundedIntegerOption(
      "branching_print_level",
      "Output level for braching code in Couenne",
      -2, J_LAST_LEVEL-1, J_WARNING,
      "");
    roptions->AddBoundedIntegerOption(
      "boundtightening_print_level",
      "Output level for bound tightening code in Couenne",
      -2, J_LAST_LEVEL-1, J_WARNING,
      "");
    roptions->AddBoundedIntegerOption(
      "convexifying_print_level",
      "Output level for convexifying code in Couenne",
      -2, J_LAST_LEVEL-1, J_WARNING,
      "");
    roptions->AddBoundedIntegerOption(
      "problem_print_level",
      "Output level for problem manipulation code in Couenne",
      -2, J_LAST_LEVEL-1, J_WARNING,
      "");


    // copied from BonminSetup::registerMilpCutGenerators(), in
    // BonBonminSetup.cpp

    struct cutOption_ {

      const char *cgname;
      int         defaultFreq;

    } cutOption [] = {
      {(const char *) "Gomory_cuts",             0},
      {(const char *) "probing_cuts",            0},
      {(const char *) "cover_cuts",              0},
      {(const char *) "mir_cuts",                0},
      {(const char *) "2mir_cuts",               0},
      {(const char *) "flow_covers_cuts",        0},
      {(const char *) "lift_and_project_cuts",   0},
      {(const char *) "reduce_split_cuts",       0},
      {(const char *) "clique_cuts",             0},
      {NULL, 0}};

    for (int i=0; cutOption [i].cgname; i++) {

      char descr [150];

      sprintf (descr, "Frequency k (in terms of nodes) for generating %s cuts in branch-and-cut.",
	       cutOption [i].cgname);

      roptions -> AddLowerBoundedIntegerOption 
	(cutOption [i].cgname,
	 descr,
	 -100, cutOption [i].defaultFreq,
	 "If k > 0, cuts are generated every k nodes, "
	 "if -99 < k < 0 cuts are generated every -k nodes but "
	 "Cbc may decide to stop generating cuts, if not enough are generated at the root node, "
	 "if k=-99 generate cuts only at the root node, if k=0 or 100 do not generate cuts.");

      roptions->setOptionExtraInfo (cutOption [i].cgname, 5);
    }
  }



  /** Add milp cut generators according to options.*/
  void CouenneSetup::addMilpCutGenerators () {

    enum extraInfo_ {CUTINFO_NONE, CUTINFO_MIG, CUTINFO_PROBING, CUTINFO_CLIQUE};

    struct cutInfo {

      const char      *optname;
      CglCutGenerator *cglptr;
      const char      *cglId;
      enum extraInfo_  extraInfo;

    } cutList [] = {
      {(const char*)"Gomory_cuts",new CglGomory,      (const char*)"Mixed Integer Gomory",CUTINFO_MIG},
      {(const char*)"probing_cuts",new CglProbing,        (const char*) "Probing",  CUTINFO_PROBING},
      {(const char*)"mir_cuts",new CglMixedIntegerRounding2, (const char*) "Mixed Integer Rounding", 
       CUTINFO_NONE},
      {(const char*)"2mir_cuts",    new CglTwomir,         (const char*) "2-MIR",    CUTINFO_NONE},
      {(const char*)"cover_cuts",   new CglKnapsackCover,  (const char*) "Cover",    CUTINFO_NONE},
      {(const char*)"clique_cuts",  new CglClique,         (const char*) "Clique",   CUTINFO_CLIQUE},
      {(const char*)"lift_and_project_cuts",new CglLandP,(const char*)"Lift and Project",CUTINFO_NONE},
      {(const char*)"reduce_split_cuts",new CglRedSplit,(const char*) "Reduce and Split",CUTINFO_NONE},
      {(const char*)"flow_covers_cuts",new CglFlowCover,(const char*) "Flow cover cuts", CUTINFO_NONE},
      {NULL, NULL, NULL, CUTINFO_NONE}};

    int freq;

    for (int i=0; cutList [i]. optname; i++) {

      options_ -> GetIntegerValue (std::string (cutList [i]. optname), freq, "bonmin.");

      if (!freq) {

	delete cutList [i].cglptr;
	continue;
      }

      CuttingMethod cg;
      cg.frequency = freq;
      cg.cgl       = cutList [i].cglptr;
      cg.id        = std::string (cutList [i]. cglId);
      cutGenerators_.push_back (cg);

      switch (cutList [i].extraInfo) {

      case CUTINFO_MIG: {
	CglGomory *gc = dynamic_cast <CglGomory *> (cutList [i].cglptr);

	if (!gc) break;

	gc -> setLimitAtRoot(512);
	gc -> setLimit(50);
      }
	break;

      case CUTINFO_PROBING: {
	CglProbing *pc = dynamic_cast <CglProbing *> (cutList [i].cglptr);

	if (!pc) break;

	pc->setUsingObjective(1);
	pc->setMaxPass(3);
	pc->setMaxPassRoot(3);
	// Number of unsatisfied variables to look at
	pc->setMaxProbe(10);
	pc->setMaxProbeRoot(50);
	// How far to follow the consequences
	pc->setMaxLook(10);
	pc->setMaxLookRoot(50);
	pc->setMaxLookRoot(10);
	// Only look at rows with fewer than this number of elements
	pc->setMaxElements(200);
	pc->setRowCuts(3);
      }
	break;

      case CUTINFO_CLIQUE: {
	CglClique *clique = dynamic_cast <CglClique *> (cutList [i].cglptr);

	if (!clique) break;

	clique -> setStarCliqueReport(false);
	clique -> setRowCliqueReport(false);
	clique -> setMinViolation(0.1);
      }
	break;

	//case CUTINFO_NONE:
      default:
	break;
      }
    }
  }
}

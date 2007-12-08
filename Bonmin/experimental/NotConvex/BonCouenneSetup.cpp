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
#include "BonCouenneInterface.hpp"

#include "CouenneObject.hpp"
#include "CouenneChooseVariable.hpp"
#include "CouenneChooseStrong.hpp"
#include "CouenneSolverInterface.hpp"
#include "CouenneCutGenerator.hpp"
#include "BonAuxInfos.hpp"
#include "BonCbcNode.hpp"

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
    //    if (CouennePtr_)
    //      delete CouennePtr_;
  }
  
  void CouenneSetup::InitializeCouenne(char **& argv){
    /* Get the basic options. */
    readOptionsFile();
    
    /** Change default value for failure behavior so that code doesn't crash when Ipopt does not solve a sub-problem.*/
    options_->SetStringValue("nlp_failure_behavior","fathom","bonmin.");

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
    int ivalue, num_points;
    options()->GetEnumValue("convexification_type", ivalue,"bonmin.");
    options()->GetIntegerValue("convexification_points",num_points,"bonmin.");
    
    CouenneCutGenerator * couenneCg = 
      new CouenneCutGenerator(ci, this, aslfg_->asl, journalist());

    CouenneProblem * couenneProb = couenneCg -> Problem();

    Bonmin::BabInfo * extraStuff = new Bonmin::BabInfo(0);

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

    // In case there are no discrete variables, we can set the optimal
    // value from the initialSolve as cutoff

    // TODO: In case there are integer variables, check if all feasible

    // pbelotti: don't set it here, it gives numerical and, with
    // Anstreicher's maximization problems, sign (!) troubles

    /*
    if (ci -> getNumIntegers() == 0) {
      // add small portion in abs. value (assume minimization problem)
      doubleParam_[Cutoff] = ci->getObjValue() + 1e-4 * fabs (ci -> getObjValue ());
    }
    */

    if(extraStuff->infeasibleNode()){
      std::cout<<"Initial linear relaxation constructed by Couenne is infeasible, quit"<<std::endl;
      return;
    }

    continuousSolver_->findIntegersAndSOS(false);

    // add CouenneObjects for branching /////////////////////////////////////////////

    {
      int nAuxs = 0, nobj = 0,
	  nVars = couenneProb -> nVars ();

      // Count # auxiliary variables

      for (int i = 0; i < nVars; i++) { // for each aux variable

	exprVar *var = couenneProb -> Var (i);

        // if this variable is associated with a nonlinear function
	if ((var -> Type () == AUX) && 
	    (var -> Image () -> Linearity () > LINEAR)) 
	  nAuxs++;
      }

      OsiObject ** objects = new OsiObject* [nAuxs];

      for (int i = 0; i < nVars; i++) { // for each aux variable

	exprVar *var = couenneProb -> Var (i);

        // if this variable is associated with a nonlinear function
	if ((var -> Type  () == AUX) && 
	    (var -> Image () -> Linearity () > LINEAR) &&
	    (var -> Multiplicity () > 0)) {

	  exprAux *aux = dynamic_cast <exprAux *> (var);

	  /*printf ("creating CouenneObject for ");
	  aux ->             print (std::cout); printf (" := ");
	  aux -> Image () -> print (std::cout); printf ("\n");*/

	  // then we may have to branch on it
	  objects [nobj] = new CouenneObject (aux, this, journalist());
	  objects [nobj++] -> setPriority (1);
	}
      }

      continuousSolver_ -> addObjects (nobj, objects);

      for(int i = 0 ; i < nobj ; i++){
       	delete objects[i];
      }

      delete [] objects;
    }

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
      dynamic_cast <CouenneSolverInterface *> (continuousSolver_) -> setCutGenPtr (couenneCg);
    }

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
    options_->GetEnumValue("varselect_stra",varSelection,"bonmin.");
    switch (varSelection) {

      // strong branching choosevariable

    case OsiTMINLPInterface::OSI_STRONG: {
	CouenneChooseStrong * chooseVariable = new CouenneChooseStrong(*this, couenneProb);
	chooseVariable->setTrustStrongForSolution(false);
	chooseVariable->setTrustStrongForBound(false);
	chooseVariable->setOnlyPseudoWhenTrusted(true);
	branchingMethod_ = chooseVariable;
	break;
    }

    case OsiTMINLPInterface::OSI_SIMPLE: // default choice
    default:
      branchingMethod_ = new CouenneChooseVariable (continuousSolver_, couenneProb);
      break;
    }

    if(intParam_[NumCutPasses] < 2)
    intParam_[NumCutPasses] = 2;
}
 
void CouenneSetup::registerOptions(){
  registerAllOptions(roptions());
}


void
  CouenneSetup::registerAllOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
    BabSetupBase::registerAllOptions(roptions);
    BonCbcFullNodeInfo::registerOptions(roptions);
    CouenneCutGenerator::registerOptions (roptions);

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
  }
  //  OsiTMINLPInterface * BonminAmplSetup::createOsiInterface{
  //}  
}

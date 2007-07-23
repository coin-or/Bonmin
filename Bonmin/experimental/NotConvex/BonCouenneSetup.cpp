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
  }
  
  void CouenneSetup::InitializeBonmin(char **& argv){
    /* Get the basic options. */
    readOptionsFile();
    
    
    /** Change default value for failure behavior so that code doesn't crash when Ipopt does not solve a sub-problem.*/
    options_->SetStringValue("nlp_failure_behavior","fathom","bonmin.");

    gatherParametersValues(options_);
    
    continuousSolver_ = new OsiClpSolverInterface;
    CouenneInterface * ci = new CouenneInterface;
    nonlinearSolver_ = ci;
    /* Read the model in various places. */
    ci->readAmplNlFile(argv,roptions(),options(),journalist());
    aslfg_ = new SmartAsl;
    aslfg_->asl = readASLfg (argv);
    
    
    /* Initialize Couenne cut generator.*/
    int ivalue, num_points;
    options()->GetEnumValue("convexification_type", ivalue,"bonmin.");
    enum conv_type convtype((enum conv_type) ivalue);
    options()->GetIntegerValue("convexification_points",num_points,"bonmin.");
    
    CouenneCutGenerator * couenneCg = new CouenneCutGenerator(ci, aslfg_->asl, true, convtype,num_points);
    CouenneProblem * couenneProb = couenneCg -> Problem();

    Bonmin::BabInfo * extraStuff = new Bonmin::BabInfo(0);
    
    // as per instructions by John Forrest, to get changed bounds
    extraStuff -> setExtraCharacteristics (extraStuff -> extraCharacteristics () | 2);
    
    continuousSolver_ -> setAuxiliaryInfo (extraStuff);
    delete extraStuff;
    extraStuff = dynamic_cast<Bonmin::BabInfo *>(continuousSolver_ -> getAuxiliaryInfo());
    /* Setup log level*/
    int lpLogLevel;
    options()->GetIntegerValue("lp_log_level",lpLogLevel,"bonmin.");
    continuousSolver_->messageHandler()->setLogLevel(lpLogLevel);
    ci->extractLinearRelaxation(*continuousSolver_, *couenneCg);
    
    if(extraStuff->infeasibleNode()){
      std::cout<<"Initial linear relaxation constructed by Couenne is infeasible, quit"<<std::endl;
      return;
    }
 
    
    continuousSolver_->findIntegersAndSOS(false);

    // [Pietro: Never used]
    //    int numberIntegerObjects = continuousSolver_->numberObjects() > 0; 

    {
      int 
	nOrig = couenneProb -> nOrig (),
	nVars = couenneProb -> nVars (),
	nAuxs = nVars - nOrig,
	nobj  = 0;

      OsiObject ** objects = new OsiObject* [nAuxs];
      
      for (int i = 0; i < nVars; i++) { // for each aux variable

	exprVar *var = couenneProb -> Var (i);

        // if this variable is associated with a nonlinear function
	if ((var -> Type () == AUX) && (var -> Image () -> Linearity () > LINEAR)) {

	  exprAux *aux = dynamic_cast <exprAux *> (var);

	  /*printf ("creating CouenneObject for ");

	  aux ->             print (std::cout); printf (" := ");
	  aux -> Image () -> print (std::cout); printf ("\n");*/

	  // then we may have to branch on it
	  objects [nobj] = new CouenneObject (aux);
	  objects [nobj++] -> setPriority (1);
	}
      }

      continuousSolver_ -> addObjects (nobj, objects);
      for(int i = 0 ; i < nobj ; i++){
       	delete objects[i];
      }
      delete [] objects;
    }
    
    //Setup Convexifier generators


    // [Pietro: never used]
    //    int numGen = 0;
    int freq;
    options()->GetIntegerValue("convexification_cuts",freq,"couenne.");
    if (freq != 0) {
      CuttingMethod cg;
      cg.frequency = freq;
      cg.cgl = couenneCg;
      cg.id = "Couenne convexifier cuts";
      cutGenerators().push_back(cg);
    }

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
      nlpHeuristic->setMaxNlpInf(1e10);
      nlpHeuristic->setNumberSolvePerLevel(numSolve);
      heuristics_.push_back(nlpHeuristic);
    }
    
    branchingMethod_ = new CouenneChooseVariable(continuousSolver_, 
                                  const_cast<CouenneProblem *> (couenneProb));

    if(intParam_[NumCutPasses] < 2)
    intParam_[NumCutPasses] = 2;
}
 
void CouenneSetup::registerOptions(){
  registerAllOptions(roptions());
}


void
  CouenneSetup::registerAllOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions){
    BabSetupBase::registerAllOptions(roptions);
    BonCbcFullNodeInfo::registerOptions(roptions);

    roptions->SetRegisteringCategory("Couenne options");
    
    roptions->AddLowerBoundedIntegerOption("convexification_cuts",
                                           "Specify the frequency (in terms of nodes) at which couenne ecp cuts are generated.",
                                           0,1,
                                           "A frequency of 0 amounts to to never solve the NLP relaxation.");
    
    roptions->AddStringOption2("local_optimization_heuristic",
                               "Do we search for local solutions of NLP's",
                               "yes",
                               "no","",
                               "yes","");
    
    roptions->AddLowerBoundedIntegerOption("log_num_local_optimization_per_level",
                               "Specify the logarithm of the number of local optimization to perform on average for each level of given depth of the tree.",
                               -1,-1,"Solve as many nlp's at the nodes for each level of the tree. Nodes are randomly selected. If for a"
                                           "given level there are less nodes than this number nlp are solved for every nodes."
                                           "For example if parameter is 8, nlp's are solved for all node untill level 8, then for half the node at level 9, 1/4 at level 10...."
                                           "Value -1 specify to perform at all nodes."
                                           );
    
    
    roptions->AddStringOption3("convexification_type",
                               "Deterimnes in which point the linear over/under-estimator are generated",
                               "current-point-only",
                               "current-point-only","Only at current optimum of relaxation",
                               "uniform-grid","Points chosen in a unform grid between the bounds of the problem",
                               "around-current-point","At points around current optimum of relaxation");
    
    roptions->AddLowerBoundedIntegerOption("convexification_points",
                                           "Specify the number of points at which to convexify when convexification type"
                                           "is uniform-grid or arround-current-point.",
                                           0,1,
                                           "");
    
    
  }
  //  OsiTMINLPInterface * BonminAmplSetup::createOsiInterface{
  //}
  
}

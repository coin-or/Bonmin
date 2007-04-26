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


#include "asl.h"
#include "getstub.h"

namespace Bonmin{
  
  SmartAsl::~SmartAsl(){
    if(asl != NULL)
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
    CouenneCutGenerator * couenneCg = new CouenneCutGenerator(ci, aslfg_->asl, true, CURRENT_ONLY,1);
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
    int numberIntegerObjects = continuousSolver_->numberObjects() > 0;
    {
      int numAuxs = couenneProb->nAuxs();
      OsiObject ** objects = new OsiObject*[numAuxs];
      int nobj = 0;
      
      for (int i = 0 ; i < numAuxs; i++) // for each aux variable
        
        // if this variable is associated with a nonlinear function
        //	if (couenneProb -> Aux (i) -> Image () -> Linearity () > LINEAR) 
      {
        /*
         printf ("creating CouenneObject for ");
         
         couenneProb -> Aux (i) ->             print (std::cout); printf (" := ");
         couenneProb -> Aux (i) -> Image () -> print (std::cout); printf ("\n");
         */
        // then we may have to branch on it
        objects [nobj] = new CouenneObject (couenneProb -> Aux (i));
        objects [nobj++] -> setPriority (1);
      }
      
      continuousSolver_ -> addObjects (nobj, objects);
      for(int i = 0 ; i < nobj ; i++){
       	delete objects[i];
      }
      delete objects;
    }
    
    //Setup Convexifier generators
    
    
    
    
    int numGen = 0;
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
    options()->GetEnumValue("nlp_local_solutions", doNlpHeurisitic, "couenne.");
    if(doNlpHeurisitic)
    {
      NlpSolveHeuristic * nlpHeuristic = new NlpSolveHeuristic;
      nlpHeuristic->setNlp(*ci,false);
      nlpHeuristic->setCouenneProblem(couenneProb);
      nlpHeuristic->setMaxNlpInf(1e10);
      heuristics_.push_back(nlpHeuristic);
    }
    
    branchingMethod_ = new CouenneChooseVariable(continuousSolver_, 
                                  const_cast<CouenneProblem *> (couenneProb));

    
}
 
void CouenneSetup::registerOptions(){
  registerAllOptions(roptions());
}


void
  CouenneSetup::registerAllOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions){
    BabSetupBase::registerAllOptions(roptions);
    roptions->SetRegisteringCategory("Couenne options");
    
    roptions->AddLowerBoundedIntegerOption("convexification_cuts",
                                           "Specify the frequency (in terms of nodes) at which couenne ecp cuts are generated.",
                                           0,1,
                                           "A frequency of 0 amounts to to never solve the NLP relaxation.");
    
    roptions->AddStringOption2("nlp_local_solutions",
                               "Do we search for local solutions of NLP's",
                               "yes",
                               "no","",
                               "yes","");
  }
  //  OsiTMINLPInterface * BonminAmplSetup::createOsiInterface{
  //}
  
}
 
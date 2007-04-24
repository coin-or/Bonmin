// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/18/2007

#include "BonCouenneSetup.hpp"
#include "BonBasicSetup.hpp"
#include "BonBonminSetup.hpp"
#include "BonNlpHeuristic.hpp"
#include "BonCouenneInterface.hpp"

#include "CouenneObject.hpp"
#include "CouenneChooseVariable.hpp"

namespace Bonmin{
  void CouenneSetup::InitializeBonmin(char **& argv){
    /* Get the basic options. */
    defaultBasicOptions();
    
    /* Read the model. */
    CouenneInterface * ci = new CouenneInterface;
    nonlinearSolver_ = ci;
    ci->readAmplNlFile(argv,journalist(),options(),roptions());
    linearSolver_ = new OsiClpSolverInterface;
    
    /* Setup log level*/
    int lpLogLevel;
    options()->GetIntegerValue("lp_log_level",lpLogLevel,"bonmin.");
    linearSolver_->messageHandler()->setLogLevel(lpLogLevel);
    ci->extractLinearRelaxation(*linearSolver_);
  
    
    
    OsiBabSolver * extraStuff = new OsiBabSolver(0);

    // as per instructions by John Forrest, to get changed bounds
    extraStuff -> setExtraCharacteristics (extraStuff -> extraCharacteristics () | 2);
    
    linearSolver_ -> setAuxiliaryInfo (extraStuff);
    delete extraStuff;
    
    linearSolver_->findIntegersAndSOS(false);
    int numberIntegerObjects = linearSolver_->numberObjects() > 0;
    {
      const CouenneProblem * couenneProb = ci->couenneProb();
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
      
      linearSolver_ -> addObjects (nobj, objects);
      //       for(int i = 0 ; i < nobj ; i++){
      // 	delete objects[i];
      //       }
      //       delete objects;
    }
    
    //Setup Convexifier generators
    
    
    
    
    int numGen = 0;
    int freq;
    options()->GetIntegerValue("convexification_cuts",freq,"couenne.");
    if (freq != 0) {
      CuttingMethod cg;
      cg.frequency = freq;
      CouenneCutGenerator *couenneGen = ci -> couenneCg ();
      cg.cgl = couenneGen;
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
      nlpHeuristic->setMaxNlpInf(1e10);
      heuristics_.push_back(nlpHeuristic);
    }
    
    branchingMethod_ = new CouenneChooseVariable(linearSolver_, 
                                  const_cast<CouenneProblem *> (ci -> couenneProb ()));

    
}
  
void
  CouenneSetup::registerAllOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions){
    BonminSetup::registerAllOptions(roptions);
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

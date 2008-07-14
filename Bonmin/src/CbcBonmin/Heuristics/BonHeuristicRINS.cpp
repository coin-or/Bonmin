// (C) Copyright CNRS and International Business Machines Corporation
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, LIF Université de la Méditérannée-CNRS
// Joao Goncalves, International Business Machines Corporation
//
// Date : 06/18/2008

#include "BonHeuristicRINS.hpp"
#include "CbcModel.hpp"
#include "OsiBranchingObject.hpp"

namespace Bonmin {

  /** Default constructor*/
  HeuristicRINS::HeuristicRINS():
    LocalSolverBasedHeuristic(){
  }
  /** Constructor with setup.*/
  HeuristicRINS::HeuristicRINS(BabSetupBase * setup):
    LocalSolverBasedHeuristic(setup){
  }

  /** Copy constructor.*/
  HeuristicRINS::HeuristicRINS
             (const HeuristicRINS &other):
    LocalSolverBasedHeuristic(other){
  }

  HeuristicRINS::~HeuristicRINS(){
  }

  /** Runs heuristic*/
  int
  HeuristicRINS::solution(double & objectiveValue,
			  double * newSolution)
  {
    if(model_->getNodeCount() || model_->getCurrentPassNumber() > 1) return 0;

    const double * bestSolution = model_->bestSolution();
    if (!bestSolution)
      return 0; // No solution found yet

    OsiTMINLPInterface * nlp = dynamic_cast<OsiTMINLPInterface *>
                               (setup_->nonlinearSolver()->clone());


    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();

    const double * currentSolution = nlp->getColSolution();

    double primalTolerance;
    nlp->getDblParam(OsiPrimalTolerance,primalTolerance);

    int nFix=0;
    for (int i=0; i<numberIntegers; i++) {
      int iColumn=integerVariable[i];
      const OsiObject * object = model_->object(i);
      // get original bounds
      double originalLower;
      double originalUpper;
      getIntegerInformation(object, originalLower, originalUpper); 
      double valueInt=bestSolution[iColumn];
      if (valueInt<originalLower) {
	valueInt=originalLower;
      } else if (valueInt>originalUpper) {
	valueInt=originalUpper;
      }
      if (fabs(currentSolution[iColumn]-valueInt)<10.0*primalTolerance) {
	double nearest=floor(valueInt+0.5);
	nlp->setColLower(iColumn,nearest);
	nlp->setColUpper(iColumn,nearest);
	nFix++;
      }
    }

    int r_val = 0;
    if(nFix > numberIntegers/5) {
      r_val = doLocalSearch(nlp, newSolution, objectiveValue, model_->getCutoff());
    }

    delete nlp;

    return r_val;
  }

  void
  HeuristicRINS::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
   roptions->SetRegisteringCategory("Local search based heuristics", RegisteredOptions::BonminCategory);
   roptions->AddStringOption2(
     "heuristic_RINS",
     "if yes runs the RINS heuristic",
     "no",
     "no", "don't run it",
     "yes", "runs the heuristic",
     "");
  }

   /** Initiaize using passed options.*/
   void 
   HeuristicRINS::Initialize(Ipopt::SmartPtr<Bonmin::OptionsList> options){
   }
}/* ends bonmin namespace*/

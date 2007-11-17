// (C) Copyright International Business Machines Corporation 2007 
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/09/2007

#include "BonNlpHeuristic.hpp"
#include "BonCouenneInterface.hpp"
#include "CouenneObject.hpp"
#include "CouenneProblem.hpp"
#include "CbcBranchActual.hpp"
#include "BonAuxInfos.hpp"

namespace Bonmin{
  NlpSolveHeuristic::NlpSolveHeuristic():
    CbcHeuristic(),
    nlp_(NULL),
    hasCloned_(false),
    maxNlpInf_(maxNlpInf_0),
    numberSolvePerLevel_(-1),
    couenne_(NULL){
  }
  
  NlpSolveHeuristic::NlpSolveHeuristic(CbcModel & model, OsiSolverInterface &nlp, bool cloneNlp, CouenneProblem * couenne):
  CbcHeuristic(model), nlp_(&nlp), hasCloned_(cloneNlp),maxNlpInf_(maxNlpInf_0),
  numberSolvePerLevel_(-1),
  couenne_(couenne){
    if(cloneNlp)
      nlp_ = nlp.clone();
  }
  
  NlpSolveHeuristic::NlpSolveHeuristic(const NlpSolveHeuristic & other):
  CbcHeuristic(other), nlp_(other.nlp_), 
  hasCloned_(other.hasCloned_),
  maxNlpInf_(other.maxNlpInf_),
  numberSolvePerLevel_(other.numberSolvePerLevel_),
  couenne_(other.couenne_){
    if(hasCloned_ && nlp_ != NULL)
      nlp_ = other.nlp_->clone();
  }
  
  CbcHeuristic * 
  NlpSolveHeuristic::clone() const{
    return new NlpSolveHeuristic(*this);
  }
  
  NlpSolveHeuristic &
  NlpSolveHeuristic::operator=(const NlpSolveHeuristic & rhs){
    if(this != &rhs){
      CbcHeuristic::operator=(rhs);
      if(hasCloned_ && nlp_)
        delete nlp_;
      
      hasCloned_ = rhs.hasCloned_;
      if(nlp_ != NULL){
        if(hasCloned_)
          nlp_ = rhs.nlp_->clone();
        else
          nlp_ = rhs.nlp_;
      }
    }
    maxNlpInf_ = rhs.maxNlpInf_;
    numberSolvePerLevel_ = rhs.numberSolvePerLevel_;
    couenne_ = rhs.couenne_;
    return *this;
  }
  
  NlpSolveHeuristic::~NlpSolveHeuristic(){
    if(hasCloned_)
      delete nlp_;
    nlp_ = NULL;
  }
  
  void
  NlpSolveHeuristic::setNlp(OsiSolverInterface &nlp, bool cloneNlp){
    if(hasCloned_ && nlp_ != NULL)
      delete nlp_;
    hasCloned_ = cloneNlp;
    if(cloneNlp)
      nlp_ = nlp.clone();
    else
      nlp_ = &nlp;
  }
  
  void
  NlpSolveHeuristic::setCouenneProblem(CouenneProblem * couenne){
    couenne_ = couenne;}
  int
  NlpSolveHeuristic::solution( double & objectiveValue, double * newSolution){
    OsiSolverInterface * solver = model_->solver();

    OsiAuxInfo * auxInfo = solver->getAuxiliaryInfo();
    BabInfo * babInfo = dynamic_cast<BabInfo *> (auxInfo);
    if(babInfo){
    babInfo->setHasNlpSolution(false);
    if(babInfo && babInfo->infeasibleNode()){
      return 0;
    }
    }

    // if too deep in the BB tree, only run NLP heuristic if
    // feasibility is low
    bool too_deep = false;

    // check depth
    if(numberSolvePerLevel_ > -1){
      if(numberSolvePerLevel_ == 0) 
	return 0;
      double num=CoinDrand48();
      const int depth = (model_->currentNode()) ? model_->currentNode()->depth() : 0;
      double beta=pow(2.,numberSolvePerLevel_);
      if(num> beta*pow(2.,-depth))
	too_deep = true;
    }

    double * lower = CoinCopyOfArray(solver->getColLower(), nlp_->getNumCols());
    double * upper = CoinCopyOfArray(solver->getColUpper(), nlp_->getNumCols());
    const double * solution = solver->getColSolution();
    OsiBranchingInformation info(solver, true);
    const int & numberObjects = model_->numberObjects();
    OsiObject ** objects = model_->objects();
    double maxInfeasibility = 0;
    for(int i = 0 ; i < numberObjects ; i++){
      CouenneObject * couObj = dynamic_cast <CouenneObject *> (objects [i]);
      if (couObj) {
	if (too_deep) { // only test infeasibility if BB level is high
	  int dummy;
	  double infeas;
	  maxInfeasibility = max ( maxInfeasibility, infeas = couObj->infeasibility(&info, dummy));
	  if(maxInfeasibility > maxNlpInf_){
	    delete [] lower;
	    delete [] upper;
	    return 0;
	  }
	}
      }
      else{
        OsiSimpleInteger * intObj = dynamic_cast<OsiSimpleInteger *>(objects[i]);
        if(intObj){
          const int & i = intObj->columnNumber();
          // Round the variable in the solver
          double value = solution[i];
          if(value < lower[i])
            value = lower[i];
          else if(value > upper[i])
            value = upper[i];
          value = floor(value + 0.5);
          lower[i] = upper[i] = value;
        }
        else{
           throw CoinError("Bonmin::NlpSolveHeuristic","solution",
                           "Unknown object.");
        }
      }
    }


    // Now set column bounds and solve the NLP with starting point
    double * saveColLower = CoinCopyOfArray(nlp_->getColLower(), nlp_->getNumCols());
    double * saveColUpper = CoinCopyOfArray(nlp_->getColUpper(), nlp_->getNumCols());
    nlp_->setColLower(lower);
    nlp_->setColUpper(upper);
    nlp_->setColSolution(solution);
    nlp_->initialSolve();
    double obj = (nlp_->isProvenOptimal()) ? nlp_->getObjValue(): DBL_MAX;

    bool foundSolution = obj < objectiveValue;
//    if(foundSolution)//Better solution found update
    {
  //    newSolution = new double [solver->getNumCols()];
      CoinCopyN(nlp_->getColSolution(), nlp_->getNumCols(), newSolution);

      //Get correct values for all auxiliary variables
      CouenneInterface * couenne = dynamic_cast<CouenneInterface *>
        (nlp_);
      if(couenne){
       couenne_->getAuxs(newSolution);
    }
      if(babInfo){
      babInfo->setNlpSolution(newSolution,model_->solver()->getNumCols(), obj);
        babInfo->setHasNlpSolution(true);
      }
      objectiveValue = obj;
  }
  nlp_->setColLower(saveColLower);
  nlp_->setColUpper(saveColUpper);
  delete [] lower;
  delete [] upper;
  delete [] saveColLower;
  delete [] saveColUpper;
  return foundSolution;
  }
}/** Ends namespace Bonmin.*/

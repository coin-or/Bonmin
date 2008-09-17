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
#include "CbcCutGenerator.hpp"
#include "CbcBranchActual.hpp"
#include "BonAuxInfos.hpp"
#include "CoinHelperFunctions.hpp"

#include "CouenneCutGenerator.hpp"
#include "CouenneProblem.hpp"

namespace Bonmin{
  NlpSolveHeuristic::NlpSolveHeuristic():
    CbcHeuristic(),
    nlp_(NULL),
    hasCloned_(false),
    maxNlpInf_(maxNlpInf_0),
    numberSolvePerLevel_(-1),
    couenne_(NULL){
    setHeuristicName("NlpSolveHeuristic");
  }
  
  NlpSolveHeuristic::NlpSolveHeuristic(CbcModel & model, OsiSolverInterface &nlp, bool cloneNlp, CouenneProblem * couenne):
  CbcHeuristic(model), nlp_(&nlp), hasCloned_(cloneNlp),maxNlpInf_(maxNlpInf_0),
  numberSolvePerLevel_(-1),
  couenne_(couenne){
    setHeuristicName("NlpSolveHeuristic");
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
      if(babInfo->infeasibleNode()){
	return 0;
      }
    }

    // if too deep in the BB tree, only run NLP heuristic if
    // feasibility is low
    bool too_deep = false;

    // check depth
    if (numberSolvePerLevel_ > -1){
      if (numberSolvePerLevel_ == 0) 
	return 0;

      const int depth = (model_ -> currentNode ()) ? model_ -> currentNode () -> depth () : 0;

      //if (CoinDrand48 () > pow (2., numberSolvePerLevel_ - depth))
      if (CoinDrand48 () > 1. / CoinMax 
	  (1., (double) ((depth - numberSolvePerLevel_) * (depth - numberSolvePerLevel_))))
	too_deep = true;
    }

    if (too_deep)
      return 0;

    double *lower = new double [couenne_ -> nVars ()];
    double *upper = new double [couenne_ -> nVars ()];

    CoinFillN (lower, couenne_ -> nVars (), -COUENNE_INFINITY);
    CoinFillN (upper, couenne_ -> nVars (),  COUENNE_INFINITY);

    CoinCopyN (solver->getColLower(), nlp_ -> getNumCols (), lower);
    CoinCopyN (solver->getColUpper(), nlp_ -> getNumCols (), upper);

    /*printf ("-- int candidate, before: ");
    for (int i=0; i<couenne_ -> nOrig (); i++) 
      printf ("[%g %g] ", lower [i], upper [i]);
      printf ("\n");*/

    const double * solution = solver->getColSolution();
    OsiBranchingInformation info (solver, true);
    const int & numberObjects = model_->numberObjects();
    OsiObject ** objects = model_->objects();
    double maxInfeasibility = 0;

    bool haveRoundedIntVars = false;

    for(int i = 0 ; i < numberObjects ; i++){
      CouenneObject * couObj = dynamic_cast <CouenneObject *> (objects [i]);
      if (couObj) {
	if (too_deep) { // only test infeasibility if BB level is high
	  int dummy;
	  double infeas;
	  maxInfeasibility = std::max ( maxInfeasibility, infeas = couObj->infeasibility(&info, dummy));
	  if(maxInfeasibility > maxNlpInf_){
	    delete [] lower;
	    delete [] upper;
	    return 0;
	  }
	}
      } else {

        OsiSimpleInteger * intObj = dynamic_cast<OsiSimpleInteger *>(objects[i]);

        if (intObj) {
          const int & i = intObj -> columnNumber ();
          // Round the variable in the solver
          double value = solution [i];
          if (value < lower[i])
	    value = lower[i];
          else if (value > upper[i])
            value = upper[i];

	  double rounded = floor (value + 0.5);

	  if (fabs (value - rounded) > COUENNE_EPS) {
	    haveRoundedIntVars = true;
	    //value = rounded;
	  }

	  // fix bounds anyway, if a better candidate is not found
	  // below at least we have an integer point
          //lower[i] = upper[i] = value;
        }
        else{
           throw CoinError("Bonmin::NlpSolveHeuristic","solution",
                           "Unknown object.");
        }
      }
    }

    // if here, it means the infeasibility is not too high. Generate a
    // better integer point as there are rounded integer variables

    bool skipOnInfeasibility = false;

    double *Y = new double [couenne_ -> nVars ()];
    CoinFillN (Y, couenne_ -> nVars (), 0.);
    CoinCopyN (solution, nlp_ -> getNumCols (), Y);

    /*printf ("-- int candidate, upon call: ");
    for (int i=0; i<couenne_ -> nOrig (); i++) 
      if (couenne_ -> Var (i) -> isInteger ())
	printf ("[%g <%g> %g] ", lower [i], Y [i], upper [i]);
      else printf ("%g ", Y [i]);
      printf ("\n");*/

    if (haveRoundedIntVars) // create "good" integer candidate for Ipopt
      skipOnInfeasibility = (couenne_ -> getIntegerCandidate (solution, Y, lower, upper) < 0);

    /*printf ("-- int candidate, after: ");
    for (int i=0; i<couenne_ -> nOrig (); i++) 
      if (couenne_ -> Var (i) -> isInteger ())
	printf ("[%g <%g> %g] ", lower [i], Y [i], upper [i]);
      else printf ("%g ", Y [i]);
      printf ("\n");*/

    bool foundSolution = false;

    if (haveRoundedIntVars && skipOnInfeasibility) 
      // no integer initial point could be found, make up some random rounding

      for (int i = couenne_ -> nOrigVars (); i--;) 

	if (couenne_ -> Var (i) -> isDefinedInteger ())
	  lower [i] = upper [i] = Y [i] = 
	    (CoinDrand48 () < 0.5) ? 
	    floor (Y [i] + COUENNE_EPS) : 
	    ceil  (Y [i] - COUENNE_EPS);

	else if (lower [i] > upper [i]) { 

	  // sanity check (should avoid problems in ex1263 with
	  // couenne.opt.obbt)

	  double swap = lower [i];
	  lower [i] = upper [i];
	  upper [i] = swap;
	}


    {
      //	printf ("[%g <%g> %g] ", lower [i], Y [i], upper [i]);

      /*printf ("int candidate: ");
	for (int i=0; i<couenne_ -> nOrig (); i++) 
	if (couenne_ -> Var (i) -> isInteger ())
	printf ("[%g <%g> %g] ", lower [i], Y [i], upper [i]);
	else printf ("%g ", Y [i]);
	printf ("\n");*/

      // Now set column bounds and solve the NLP with starting point
      double * saveColLower = CoinCopyOfArray (nlp_ -> getColLower (), nlp_ -> getNumCols ());
      double * saveColUpper = CoinCopyOfArray (nlp_ -> getColUpper (), nlp_ -> getNumCols ());

      for (int i = nlp_ -> getNumCols (); i--;) {

	if (lower [i] > upper [i]) {
	  double swap = lower [i];
	  lower [i] = upper [i];
	  upper [i] = swap;
	}

	if      (Y [i] < lower [i]) Y [i] = lower [i];
	else if (Y [i] > upper [i]) Y [i] = upper [i];
      }


      nlp_ -> setColLower    (lower);
      nlp_ -> setColUpper    (upper);
      nlp_ -> setColSolution (Y);

      // apply NLP solver /////////////////////////////////
      nlp_ -> initialSolve ();

      double obj = (nlp_ -> isProvenOptimal()) ? nlp_ -> getObjValue (): COIN_DBL_MAX;

      if (nlp_ -> isProvenOptimal () &&
	  couenne_ -> checkNLP (nlp_ -> getColSolution (), obj, true) && // true for recomputing obj
	  (obj < couenne_ -> getCutOff ())) {

	// store solution in Aux info

	const int nVars = solver->getNumCols();
	double* tmpSolution = new double [nVars];
	CoinCopyN (nlp_ -> getColSolution(), nlp_ -> getNumCols(), tmpSolution);

	//Get correct values for all auxiliary variables
	CouenneInterface * couenne = dynamic_cast <CouenneInterface *> (nlp_);

	if (couenne)
	  couenne_ -> getAuxs (tmpSolution);

	if (babInfo){
	  babInfo->setNlpSolution (tmpSolution, nVars, obj);
	  babInfo->setHasNlpSolution (true);
	}

	if (obj < objectiveValue) { // found better solution?

	  const CouNumber 
	    *lb = solver -> getColLower (),
	    *ub = solver -> getColUpper ();

	  // check bounds once more after getAux. This avoids false
	  // asserts in CbcModel.cpp:8305 on integerTolerance violated
	  for (int i=0; i < nVars; i++, lb++, ub++) {

	    CouNumber &t = tmpSolution [i];
	    if      (t < *lb) t = *lb;
	    else if (t > *ub) t = *ub;
	  }

	  couenne_ -> setCutOff (obj);
	  foundSolution = true;
	  objectiveValue = obj;
	  CoinCopyN (tmpSolution, nVars, newSolution);
	}
	delete [] tmpSolution;
      }

      nlp_->setColLower (saveColLower);
      nlp_->setColUpper (saveColUpper);

      delete [] saveColLower;
      delete [] saveColUpper;
    }

    delete [] Y;

    delete [] lower;
    delete [] upper;

    return foundSolution;
  }
}/** Ends namespace Bonmin.*/

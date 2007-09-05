// Copyright (C) 2006, 2007 International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <climits>
#include "BonChooseVariable.hpp"
#include "CoinTime.hpp"
#include "IpBlas.hpp"

// This couples Cbc code into Bonmin code...
#include "CbcModel.hpp"
#include "BonGuessHeuristic.hpp"

namespace Bonmin {

BonChooseVariable::BonChooseVariable(OsiTMINLPInterface * solver) :
  OsiChooseStrong(solver),
  cbc_model_(NULL),
  only_pseudo_when_trusted_(false),
  guessHeuristic_(NULL)
{
  SmartPtr<TNLPSolver> tnlp_solver =
    static_cast<TNLPSolver *> (solver->solver());
  DBG_ASSERT(IsValid(tnlp_solver));
  jnlst_ = tnlp_solver->Jnlst();
  DBG_ASSERT(IsValid(jnlst_));
  SmartPtr<OptionsList> options = tnlp_solver->Options();

  options->GetIntegerValue("bb_log_level", bb_log_level_, "bonmin.");
  options->GetNumericValue("setup_pseudo_frac", setup_pseudo_frac_, "bonmin.");
  options->GetNumericValue("maxmin_crit_no_sol", maxmin_crit_no_sol_, "bonmin.");
  options->GetNumericValue("maxmin_crit_have_sol", maxmin_crit_have_sol_, "bonmin.");
  if (!options->GetIntegerValue("number_before_trust_list", numberBeforeTrustedList_, "bonmin.")) {
    // default is to use the same value as for numberBeforeTrusted
    options->GetIntegerValue("number_before_trust", numberBeforeTrustedList_, "bonmin.");
  }
  options->GetIntegerValue("number_strong_branch_root", numberStrongRoot_, "bonmin.");

}

BonChooseVariable::BonChooseVariable(const BonChooseVariable & rhs) :
  OsiChooseStrong(rhs),
  cbc_model_(rhs.cbc_model_),
  only_pseudo_when_trusted_(rhs.only_pseudo_when_trusted_),
  maxmin_crit_no_sol_(rhs.maxmin_crit_no_sol_),
  maxmin_crit_have_sol_(rhs.maxmin_crit_have_sol_),
  setup_pseudo_frac_(rhs.setup_pseudo_frac_),
  numberBeforeTrustedList_(rhs.numberBeforeTrustedList_),
  numberStrongRoot_(rhs.numberStrongRoot_),
  guessHeuristic_(rhs.guessHeuristic_)
{
  jnlst_ = rhs.jnlst_;
  bb_log_level_ = rhs.bb_log_level_;
  DBG_ASSERT(IsValid(jnlst_));
}

BonChooseVariable &
BonChooseVariable::operator=(const BonChooseVariable & rhs)
{
  if (this != &rhs) {
    OsiChooseStrong::operator=(rhs);
    jnlst_ = rhs.jnlst_;
    bb_log_level_ = rhs.bb_log_level_;
    cbc_model_ = rhs.cbc_model_;
    only_pseudo_when_trusted_ = rhs.only_pseudo_when_trusted_;
    maxmin_crit_no_sol_ = rhs.maxmin_crit_no_sol_;
    maxmin_crit_have_sol_ = rhs.maxmin_crit_have_sol_;
    setup_pseudo_frac_ = rhs.setup_pseudo_frac_;
    numberBeforeTrustedList_ = rhs.numberBeforeTrustedList_;
    numberStrongRoot_ = rhs.numberStrongRoot_;
    guessHeuristic_ = rhs.guessHeuristic_;
  }
  return *this;
}

OsiChooseVariable *
BonChooseVariable::clone() const
{
  return new BonChooseVariable(*this);
}

BonChooseVariable::~BonChooseVariable ()
{}

int
BonChooseVariable::setupList ( OsiBranchingInformation *info, bool initialize)
{
  if (guessHeuristic_) {
    guessHeuristic_->setChooseMethod(this);
    //guessHeuristic_ = NULL;
  }
  if (numberBeforeTrustedList_ < 0) {
    number_not_trusted_ = 1;
    return OsiChooseVariable::setupList(info, initialize);
  }
  if (initialize) {
    status_=-2;
    delete [] goodSolution_;
    bestObjectIndex_=-1;
    numberStrongDone_=0;
    numberStrongIterations_ = 0;
    numberStrongFixed_ = 0;
    goodSolution_ = NULL;
    goodObjectiveValue_ = COIN_DBL_MAX;
    number_not_trusted_=0;
  }
  else {
    throw -1;
  }
  numberOnList_=0;
  numberUnsatisfied_=0;
  int numberObjects = solver_->numberObjects();
  assert (numberObjects);
  if (numberObjects>numberObjects_) {
    //AW : How could that ever happen?  Right now, all old content is deleted!
    assert(false && "Right now, all old content is deleted!");
    // redo useful arrays
    OsiPseudoCosts::initialize(numberObjects);
  }
  double check = -COIN_DBL_MAX;
  int checkIndex=0;
  int bestPriority=COIN_INT_MAX;
  int maximumStrong= CoinMin(CoinMax(numberStrong_,numberStrongRoot_),
			     numberObjects) ;
  int putOther = numberObjects;
  int i;
  for (i=0;i<numberObjects;i++) {
    list_[i]=-1;
    useful_[i]=0.0;
  }
  // We make a second list for most fractional variables
  int* list2 = NULL;
  double* useful2 = NULL;
  double check2 = -COIN_DBL_MAX;
  int checkIndex2=0;
  int max_most_fra = setup_pseudo_frac_ > 0. ? (int)floor(setup_pseudo_frac_*(double)maximumStrong): 0;
  if (setup_pseudo_frac_ > 0.) {
    max_most_fra = CoinMax(1, max_most_fra);
  }
  if (max_most_fra) {
    list2 = new int[max_most_fra];
    useful2 = new double[max_most_fra];
    for (i=0;i<max_most_fra;i++) {
      list2[i]=-1;
      useful2[i]=0.0;
    }
  }

  OsiObject ** object = info->solver_->objects();
  // Get average pseudo costs and see if pseudo shadow prices possible
  int shadowPossible=shadowPriceMode_;
  assert(!shadowPossible);
  if (shadowPossible) {
    for ( i=0;i<numberObjects;i++) {
      if ( !object[i]->canHandleShadowPrices()) {
	shadowPossible=0;
	break;
      }
    }
    if (shadowPossible) {
      int numberRows = solver_->getNumRows();
      const double * pi = info->pi_;
      double sumPi=0.0;
      for (i=0;i<numberRows;i++) 
	sumPi += fabs(pi[i]);
      sumPi /= ((double) numberRows);
      // and scale back
      sumPi *= 0.01;
      info->defaultDual_ = sumPi; // switch on
      int numberColumns = solver_->getNumCols();
      int size = CoinMax(numberColumns,2*numberRows);
      info->usefulRegion_ = new double [size];
      CoinZeroN(info->usefulRegion_,size);
      info->indexRegion_ = new int [size];
    }
  }
  double sumUp=0.0;
  double numberUp=0.0;
  double sumDown=0.0;
  double numberDown=0.0;
  for ( i=0;i<numberObjects;i++) {
    sumUp += upTotalChange_[i];
    numberUp += upNumber_[i];
    sumDown += downTotalChange_[i];
    numberDown += downNumber_[i];
    if (bb_log_level_>4) {
      printf("%3d up %3d  %15.8e  down %3d  %15.8e\n",i, upNumber_[i], upTotalChange_[i], downNumber_[i], downTotalChange_[i]);
    }
  }
  double upMultiplier=(1.0+sumUp)/(1.0+numberUp);
  double downMultiplier=(1.0+sumDown)/(1.0+numberDown);
  if (bb_log_level_>4) {
    printf("upMultiplier = %e downMultiplier = %e\n",upMultiplier,downMultiplier);
  }
  // Say feasible
  bool feasible = true;
  for ( i=0;i<numberObjects;i++) {
    int way;
    double value = object[i]->infeasibility(info,way);
    if (value>0.0) {
      numberUnsatisfied_++;
      if (value>=1e50) {
	// infeasible
	feasible=false;
	break;
      }
      int priorityLevel = object[i]->priority();
      // Better priority? Flush choices.
      if (priorityLevel<bestPriority) {
	for (int j=maximumStrong-1;j>=0;j--) {
	  if (list_[j]>=0) {
	    int iObject = list_[j];
	    list_[j]=-1;
	    useful_[j]=0.0;
	    list_[--putOther]=iObject;
	  }
	}
	maximumStrong = CoinMin(maximumStrong,putOther);
	bestPriority = priorityLevel;
	check=-COIN_DBL_MAX;
	checkIndex=0;
	check2=-COIN_DBL_MAX;
	checkIndex2=0;
	number_not_trusted_=0;
	if (max_most_fra>0) {
	  for (int j=0;j<max_most_fra;j++) {
	    list2[j]=-1;
	    useful2[j]=0.0;
	  }
	}
      } 
      if (priorityLevel==bestPriority) {
	// Modify value
	sumUp = upTotalChange_[i]+1.0e-30;
	numberUp = upNumber_[i];
	sumDown = downTotalChange_[i]+1.0e-30;
	numberDown = downNumber_[i];
	double upEstimate = object[i]->upEstimate();
	double downEstimate = object[i]->downEstimate();
	if (shadowPossible<2) {
	  upEstimate = numberUp ? ((upEstimate*sumUp)/numberUp) : (upEstimate*upMultiplier);
	  //if (numberUp<numberBeforeTrusted_)
	  //  upEstimate *= (numberBeforeTrusted_+1.0)/(numberUp+1.0);
	  downEstimate = numberDown ? ((downEstimate*sumDown)/numberDown) : (downEstimate*downMultiplier);
	  //if (numberDown<numberBeforeTrusted_)
	  //  downEstimate *= (numberBeforeTrusted_+1.0)/(numberDown+1.0);
	} else {
	  // use shadow prices always
	}
	double value2 = -COIN_DBL_MAX;
	if (numberUp<numberBeforeTrustedList_ ||
	    numberDown<numberBeforeTrustedList_) {
	  value2 = value;
	}
	double MAXMIN_CRITERION = maxminCrit();
	value = MAXMIN_CRITERION*CoinMin(upEstimate,downEstimate) + (1.0-MAXMIN_CRITERION)*CoinMax(upEstimate,downEstimate);
	if (bb_log_level_>4) {
	  printf("%3d value = %e upEstimate = %e downEstimate = %e infeas = %e value2 = %e\n", i,value,upEstimate,downEstimate,object[i]->infeasibility(info,way),value2);
	}
	if (value>check) {
	  //add to list
	  int iObject = list_[checkIndex];
	  if (iObject>=0) {
	    assert (list_[putOther-1]<0);
	    list_[--putOther]=iObject;  // to end
	  }
	  list_[checkIndex]=i;
	  assert (checkIndex<putOther);
	  useful_[checkIndex]=value;
	  // find worst
	  check=COIN_DBL_MAX;
	  maximumStrong = CoinMin(maximumStrong,putOther);
	  for (int j=0;j<maximumStrong;j++) {
	    if (list_[j]>=0) {
	      if (useful_[j]<check) {
		check=useful_[j];
		checkIndex=j;
	      }
	    } else {
	      check=0.0;
	      checkIndex = j;
	      break;
	    }
	  }
	} else {
	  // to end
	  assert (list_[putOther-1]<0);
	  list_[--putOther]=i;
	  maximumStrong = CoinMin(maximumStrong,putOther);
	}
	if (max_most_fra > 0 && value2>check2) {
	  // add to list of integer infeasibilities
	  number_not_trusted_++;
	  list2[checkIndex2]=i;
	  useful2[checkIndex2]=value2;
	  // find worst
	  check2=COIN_DBL_MAX;
	  for (int j=0;j<max_most_fra;j++) {
	    if (list2[j]>=0) {
	      if (useful2[j]<check2) {
		check2=useful2[j];
		checkIndex2=j;
	      }
	    } else {
	      check2=0.0;
	      checkIndex2 = j;
	      break;
	    }
	  }
	}
      } else {
	// worse priority
	// to end
	assert (list_[putOther-1]<0);
	list_[--putOther]=i;
	maximumStrong = CoinMin(maximumStrong,putOther);
      }
    }
  }
#if 0
  for (int i=0; i<maximumStrong; i++) { int way;
      printf("list_[%5d] = %5d, usefull_[%5d] = %23.16e %23.16e \n", i,list_[i],i,useful_[i],object[list_[i]]->infeasibility(info,way));
  }
#endif
  // Get list
  numberOnList_=0;
  if (feasible) {
    maximumStrong = CoinMin(maximumStrong,putOther);
    for (i=0;i<maximumStrong;i++) {
      if (list_[i]>=0) {
	list_[numberOnList_]=list_[i];
	useful_[numberOnList_++]=-useful_[i];
	if (bb_log_level_>4) {
	  printf("list_[%3d] = %3d useful_[%3d] = %e\n",numberOnList_-1,list_[numberOnList_-1],numberOnList_-1,useful_[numberOnList_-1]);
	}
      }
    }
    if (numberOnList_) {
      int tmp_on_list = 0;
      if (max_most_fra > 0 && numberOnList_ >= maximumStrong) {
	// If we want to force non-trusted in the list, give them huge
	// weight here
	number_not_trusted_=0;
	for (i=0;i<max_most_fra;i++) {
	  if (list2[i]>=0) {
	    list2[number_not_trusted_] = list2[i];
	    useful2[number_not_trusted_++] = useful2[i];
	    if (bb_log_level_>4) {
	      printf("list2[%3d] = %3d useful2[%3d] = %e\n",number_not_trusted_-1,list2[number_not_trusted_-1],number_not_trusted_-1,useful2[number_not_trusted_-1]);
	    }
	  }
	}
	if (number_not_trusted_) {
	  CoinSort_2(list_,list_+numberOnList_,useful_);
	  CoinSort_2(list2,list2+number_not_trusted_,useful2);
	  int i1=0;
	  int i2=0;
	  for (i=0; i<numberObjects; i++) {
	    bool found1 = (list_[i1]==i);
	    bool found2 = (list2[i2]==i);
	    if (found1 && found2) {
	      useful_[i1] = -1e150*(1.+useful2[i2]);
	      list2[i2] = -1;
	    }
	    if (found1) i1++;
	    if (found2) i2++;
	    if (i2==max_most_fra) break;
	  }
	  for (i=0; i<number_not_trusted_; i++) {
	    if (list2[i] >= 0) {
	      list_[numberOnList_+tmp_on_list] = list2[i];
	      useful_[numberOnList_+tmp_on_list] = -1e150*(1.+useful2[i]);
	      tmp_on_list++;
	    }
	  }
	}
      }
      // Sort 
      CoinSort_2(useful_,useful_+numberOnList_+tmp_on_list,list_);
      // move others
      i = numberOnList_;
      for (;putOther<numberObjects;putOther++) 
	list_[i++]=list_[putOther];
      assert (i==numberUnsatisfied_);
      if (!CoinMax(numberStrong_,numberStrongRoot_))
	numberOnList_=0;
    }
  } else {
    // not feasible
    numberUnsatisfied_=-1;
  }
  // Get rid of any shadow prices info
  info->defaultDual_ = -1.0; // switch off
  delete [] info->usefulRegion_;
  delete [] info->indexRegion_;
  delete [] list2;
  delete [] useful2;
  int way;
  if (bb_log_level_>3) {
    //for (int i=0; i<Min(numberUnsatisfied_,numberStrong_); i++)
    for (int i=0; i<numberOnList_; i++)
      printf("list_[%5d] = %5d, usefull_[%5d] = %23.16e %23.16e \n", i,list_[i],i,useful_[i],object[list_[i]]->infeasibility(info,way));
  }
  return numberUnsatisfied_;

}

/* Choose a variable
   Returns - 
   -1 Node is infeasible
   0  Normal termination - we have a candidate
   1  All looks satisfied - no candidate
   2  We can change the bound on a variable - but we also have a strong branching candidate
   3  We can change the bound on a variable - but we have a non-strong branching candidate
   4  We can change the bound on a variable - no other candidates
   We can pick up branch from whichObject() and whichWay()
   We can pick up a forced branch (can change bound) from whichForcedObject() and whichForcedWay()
   If we have a solution then we can pick up from goodObjectiveValue() and goodSolution()
*/
int 
BonChooseVariable::chooseVariable(
  OsiSolverInterface * solver,
  OsiBranchingInformation *info,
  bool fixVariables)
{
  // We assume here that chooseVariable is called once at the very
  // beginning with fixVariables set to true.  This is then the root
  // node.
  bool isRoot = isRootNode();
  if (!isRoot) {
    numberStrongRoot_ = numberStrong_;
  }
  if (numberUnsatisfied_) {
    int numberLeft = CoinMin(numberStrongRoot_-numberStrongDone_,numberUnsatisfied_);
    int numberToDo=0;
    int * temp = (int *) useful_;
    OsiHotInfo * results = new OsiHotInfo [numberLeft];
    int returnCode=0;
    bestObjectIndex_ = -1;
    bestWhichWay_ = -1;
    firstForcedObjectIndex_ = -1;
    firstForcedWhichWay_ =-1;
    double bestTrusted=-COIN_DBL_MAX;
    for (int i=0;i<numberLeft;i++) {
      int iObject = list_[i];
      if (numberBeforeTrusted_<0||
	  (only_pseudo_when_trusted_ && number_not_trusted_>0) ||
	  !isRoot && (upNumber_[iObject]<numberBeforeTrusted_ ||
		      downNumber_[iObject]<numberBeforeTrusted_ )||
	  isRoot && (!upNumber_[iObject] && !downNumber_[iObject]) ) {
	results[numberToDo] = OsiHotInfo(solver,info,const_cast<const OsiObject **> (solver->objects()),iObject);
	temp[numberToDo++]=iObject;
      } else {
	const OsiObject * obj = solver->object(iObject);
	double upEstimate = (upTotalChange_[iObject]*obj->upEstimate())/upNumber_[iObject];
	double downEstimate = (downTotalChange_[iObject]*obj->downEstimate())/downNumber_[iObject];
	double MAXMIN_CRITERION = maxminCrit();
	double value = MAXMIN_CRITERION*CoinMin(upEstimate,downEstimate) + (1.0-MAXMIN_CRITERION)*CoinMax(upEstimate,downEstimate);
	if (value > bestTrusted) {
	  bestObjectIndex_=iObject;
	  bestWhichWay_ = upEstimate>downEstimate ? 0 : 1;
	  bestTrusted = value;
	}
      }
    }
    int numberFixed=0;
    if (numberToDo) {
      int numberDone=0;
      returnCode = doBonStrongBranching(solver,info,numberToDo,results,1,numberDone);
      if (bb_log_level_>=3) {
	const char* stat_msg[] = {"NOTDON", "FEAS", "INFEAS", "NOFINI"};
	printf("BON0001I           DownStat    DownChange     UpStat      UpChange\n");
	for (int i = 0; i<numberDone; i++) {
	  double up_change = results[i].upChange();
	  double down_change = results[i].downChange();
	  int up_status = results[i].upStatus();
	  int down_status = results[i].downStatus();
	  printf("BON0002I    %3d    %6s    %13.6e   %6s    %13.6e\n",
		 i, stat_msg[down_status+1], down_change, stat_msg[up_status+1], up_change);
	}
      }
      if (returnCode>=0&&returnCode<=2) {
	if (returnCode) {
	  if (returnCode==2)
	    numberToDo = numberDone;
	  returnCode=4;
	  if (bestObjectIndex_>=0)
	    returnCode=3;
	}
	for (int i=0;i<numberToDo;i++) {
	  int iObject = temp[i];
	  double upEstimate;
	  if (results[i].upStatus()!=1) {
	    assert (results[i].upStatus()>=0);
	    upEstimate = results[i].upChange();
	  } else {
	    // infeasible - just say expensive
	    if (info->cutoff_<1.0e50)
	      upEstimate = 2.0*(info->cutoff_-info->objectiveValue_);
	    else
	      upEstimate = 2.0*fabs(info->objectiveValue_);
	    if (firstForcedObjectIndex_ <0) {
	      // first fixed variable
	      firstForcedObjectIndex_ = iObject;
	      firstForcedWhichWay_ =0;
	    }
	    numberFixed++;
	    if (fixVariables) {
	      const OsiObject * obj = solver->object(iObject);
	      OsiBranchingObject * branch = obj->createBranch(solver,info,0);
	      branch->branch(solver);
	      delete branch;
	    }
	  }
	  double downEstimate;
	  if (results[i].downStatus()!=1) {
	    assert (results[i].downStatus()>=0);
	    downEstimate = results[i].downChange();
	  } else {
	    // infeasible - just say expensive
	    if (info->cutoff_<1.0e50)
	      downEstimate = 2.0*(info->cutoff_-info->objectiveValue_);
	    else
	      downEstimate = 2.0*fabs(info->objectiveValue_);
	    if (firstForcedObjectIndex_ <0) {
	      firstForcedObjectIndex_ = iObject;
	      firstForcedWhichWay_ =1;
	    }
	    numberFixed++;
	    if (fixVariables) {
	      const OsiObject * obj = solver->object(iObject);
	      OsiBranchingObject * branch = obj->createBranch(solver,info,1);
	      branch->branch(solver);
	      delete branch;
	    }
	  }
	  double MAXMIN_CRITERION = maxminCrit();
	  double value = MAXMIN_CRITERION*CoinMin(upEstimate,downEstimate) + (1.0-MAXMIN_CRITERION)*CoinMax(upEstimate,downEstimate);
	  if (value>bestTrusted) {
	    bestTrusted = value;
	    bestObjectIndex_ = iObject;
	    bestWhichWay_ = upEstimate>downEstimate ? 0 : 1;
	    // but override if there is a preferred way
	    const OsiObject * obj = solver->object(iObject);
	    if (obj->preferredWay()>=0&&obj->infeasibility())
	      bestWhichWay_ = obj->preferredWay();
	    if (returnCode)
	      returnCode=2;
	  }
	}
      } else if (returnCode==3) {
	// max time - just choose one
	bestObjectIndex_ = list_[0];
	bestWhichWay_ = 0;
	returnCode=0;
      }
    } else {
      bestObjectIndex_=list_[0];
    }
    delete [] results;
    if ( bestObjectIndex_ >=0 ) {
      OsiObject * obj = solver->objects()[bestObjectIndex_];
      obj->setWhichWay(	bestWhichWay_);
    }
    if (bb_log_level_>4) {
      printf("           Choosing %d\n", bestObjectIndex_);
    }
    if (numberFixed==numberUnsatisfied_&&numberFixed)
      returnCode=4;
    return returnCode;
  } else {
    return 1;
  }

}

int 
BonChooseVariable::doBonStrongBranching( OsiSolverInterface * solver, 
					 OsiBranchingInformation *info, int numberToDo,
					 OsiHotInfo * results, int returnCriterion,
					 int & numberDone)
{
  // Might be faster to extend branch() to return bounds changed
  double * saveLower = NULL;
  double * saveUpper = NULL;
  int numberColumns = solver->getNumCols();
  solver->markHotStart();
  const double * lower = info->lower_;
  const double * upper = info->upper_;
  saveLower = CoinCopyOfArray(info->lower_,numberColumns);
  saveUpper = CoinCopyOfArray(info->upper_,numberColumns);
  numberDone=0;
  int returnCode=0;
  double timeStart = CoinCpuTime();
  for (int iDo=0;iDo<numberToDo;iDo++) {
    OsiHotInfo * result = results + iDo;
    // For now just 2 way
    OsiBranchingObject * branch = result->branchingObject();
    assert (branch->numberBranches()==2);
    /*
      Try the first direction.  Each subsequent call to branch() performs the
      specified branch and advances the branch object state to the next branch
      alternative.)
    */
    OsiSolverInterface * thisSolver = solver; 
    if (branch->boundBranch()) {
      // ordinary
      branch->branch(solver);
      // maybe we should check bounds for stupidities here?
      solver->solveFromHotStart() ;
    } else {
      // adding cuts or something 
      thisSolver = solver->clone();
      branch->branch(thisSolver);
      // set hot start iterations
      int limit;
      thisSolver->getIntParam(OsiMaxNumIterationHotStart,limit);
      thisSolver->setIntParam(OsiMaxNumIteration,limit); 
      thisSolver->resolve();
    }
    // can check if we got solution
    // status is 0 finished, 1 infeasible and 2 unfinished and 3 is solution
    int status0 = result->updateInformation(thisSolver,info,this);
    // is decrease already below, once should be enough... ??????
    numberStrongDone_++;
    numberStrongIterations_ += thisSolver->getIterationCount();
    if (status0==3) {
      // new solution already saved
      if (trustStrongForSolution_) {
	info->cutoff_ = goodObjectiveValue_;
	status0=0;
      }
    }
    if (solver!=thisSolver)
      delete thisSolver;
    // Restore bounds
    for (int j=0;j<numberColumns;j++) {
      if (saveLower[j] != lower[j])
	solver->setColLower(j,saveLower[j]);
      if (saveUpper[j] != upper[j])
	solver->setColUpper(j,saveUpper[j]);
    }
    /*
      Try the next direction
    */
    thisSolver = solver; 
    if (branch->boundBranch()) {
      // ordinary
      branch->branch(solver);
      // maybe we should check bounds for stupidities here?
      solver->solveFromHotStart() ;
    } else {
      // adding cuts or something 
      thisSolver = solver->clone();
      branch->branch(thisSolver);
      // set hot start iterations
      int limit;
      thisSolver->getIntParam(OsiMaxNumIterationHotStart,limit);
      thisSolver->setIntParam(OsiMaxNumIteration,limit); 
      thisSolver->resolve();
    }
    // can check if we got solution
    // status is 0 finished, 1 infeasible and 2 unfinished and 3 is solution
    int status1 = result->updateInformation(thisSolver,info,this);
    numberStrongDone_++;
    numberStrongIterations_ += thisSolver->getIterationCount();
    if (status1==3) {
      // new solution already saved
      if (trustStrongForSolution_) {
	info->cutoff_ = goodObjectiveValue_;
	status1=0;
      }
    }
    if (solver!=thisSolver)
      delete thisSolver;
    // Restore bounds
    for (int j=0;j<numberColumns;j++) {
      if (saveLower[j] != lower[j])
	solver->setColLower(j,saveLower[j]);
      if (saveUpper[j] != upper[j])
	solver->setColUpper(j,saveUpper[j]);
    }

    /*
      End of evaluation for this candidate variable. Possibilities are:
      * Both sides below cutoff; this variable is a candidate for branching.
      * Both sides infeasible or above the objective cutoff: no further action
      here. Break from the evaluation loop and assume the node will be purged
      by the caller.
      * One side below cutoff: Install the branch (i.e., fix the variable). Possibly break
      from the evaluation loop and assume the node will be reoptimised by the
      caller.
    */
    numberDone++;
    if (status0==1&&status1==1) {
      // infeasible
      returnCode=-1;
      break; // exit loop
    } else if (status0==1||status1==1) {
      numberStrongFixed_++;
      if (!returnCriterion) {
	returnCode=1;
      } else {
	returnCode=2;
	break;
      }
    }
    bool hitMaxTime = ( CoinCpuTime()-timeStart > info->timeRemaining_);
    if (hitMaxTime) {
      returnCode=3;
      break;
    }
  }
  delete [] saveLower;
  delete [] saveUpper;
  // Delete the snapshot
  solver->unmarkHotStart();
  return returnCode;
}

bool BonChooseVariable::isRootNode() const {
  assert(cbc_model_);
  const int depth = cbc_model_->currentNode()->depth();
  return (depth == 0);
}

double
BonChooseVariable::maxminCrit() const {
  double retval = maxmin_crit_no_sol_;
  if (cbc_model_) {
    const int stateOfSearch = cbc_model_->stateOfSearch();
    const int depth = cbc_model_->currentNode()->depth();
    if (stateOfSearch>1 && depth>10) {
      retval = maxmin_crit_have_sol_;
    }
  }
  return retval;
}

// Given a candidate  fill in useful information e.g. estimates
void 
BonChooseVariable::updateInformation(const OsiBranchingInformation *info,
				     int branch, OsiHotInfo * hotInfo)
{
  int index = hotInfo->whichObject();
  assert (index<solver_->numberObjects());
  const OsiObject * object = info->solver_->object(index);
  assert (object->upEstimate()>0.0&&object->downEstimate()>0.0);
  assert (branch<2);
  if (branch) {
    //if (hotInfo->upStatus()!=1) {
    // AW: Let's update the pseudo costs only if the strong branching
    // problem was marked as "solved"
    if (hotInfo->upStatus()==0) {
      assert (hotInfo->upStatus()>=0);
      upTotalChange_[index] += hotInfo->upChange()/object->upEstimate();
      upNumber_[index]++;
    } else if (hotInfo->upStatus()==1) {
      // infeasible - just say expensive
      upNumber_[index]++;
      if (info->cutoff_<1.0e50)
	upTotalChange_[index] += 2.0*(info->cutoff_-info->objectiveValue_)/object->upEstimate();
      else
	upTotalChange_[index] += 2.0*fabs(info->objectiveValue_)/object->upEstimate();
    }
  } else {
    if (hotInfo->downStatus()==0) {
      assert (hotInfo->downStatus()>=0);
      downTotalChange_[index] += hotInfo->downChange()/object->downEstimate();
      downNumber_[index]++;
    } else if (hotInfo->upStatus()==1) {
      downNumber_[index]++;
      // infeasible - just say expensive
      if (info->cutoff_<1.0e50)
	downTotalChange_[index] += 2.0*(info->cutoff_-info->objectiveValue_)/object->downEstimate();
      else
	downTotalChange_[index] += 2.0*fabs(info->objectiveValue_)/object->downEstimate();
    }
  }  
}
// Given a branch fill in useful information e.g. estimates
void 
BonChooseVariable::updateInformation( int index, int branch, 
				      double changeInObjective, double changeInValue,
				      int status)
{
  assert (index<solver_->numberObjects());
  assert (branch<2);
  assert (changeInValue>0.0);
  assert (branch<2);
  //printf("update %3d %3d %e %e %3d\n",index, branch, changeInObjective,changeInValue,status);
  if (branch) {
    if (status!=1) {
      assert (status>=0);
      upTotalChange_[index] += changeInObjective/changeInValue;
      upNumber_[index]++;
    } else {
      // infeasible - just say expensive
      assert(cbc_model_); // Later, we need to get this information in a different way...
      upNumber_[index]++;
      double cutoff = cbc_model_->getCutoff();
      double objectiveValue = cbc_model_->getCurrentObjValue();
      if (cutoff<1.0e50)
	upTotalChange_[index] += 2.0*(cutoff-objectiveValue)/changeInValue;
      else
	upTotalChange_[index] += 2.0*fabs(objectiveValue)/changeInValue;
    }
  } else {
    if (status!=1) {
      assert (status>=0);
      downTotalChange_[index] += changeInObjective/changeInValue;
      downNumber_[index]++;
    } else {
      assert(cbc_model_);
      // infeasible - just say expensive
      downNumber_[index]++;
      double cutoff = cbc_model_->getCutoff();
      double objectiveValue = cbc_model_->getCurrentObjValue();
      if (cutoff<1.0e50)
	downTotalChange_[index] += 2.0*(cutoff-objectiveValue)/changeInValue;
      else
	downTotalChange_[index] += 2.0*fabs(objectiveValue)/changeInValue;
    }
  }  
}

}/* Ends Bonmin's namespace.*/

/*
 * Name:    CouenneChooseStrong.cpp
 * Authors: Andreas Waechter, IBM Corp.
 * Purpose: Strong branching objects for Couenne
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneChooseStrong.hpp"

CouenneChooseStrong::CouenneChooseStrong(Bonmin::BabSetupBase &b,
					 CouenneProblem* si):
  OsiChooseStrong(const_cast<Bonmin::OsiTMINLPInterface *>(b.nonlinearSolver())),
  cbc_model_(NULL),
  only_pseudo_when_trusted_(false),
  problem_(si)
{
  jnlst_ = b.journalist();
  SmartPtr<OptionsList> options = b.options();

  options->GetIntegerValue("bb_log_level", bb_log_level_, "bonmin.");
  options->GetNumericValue("setup_pseudo_frac", setup_pseudo_frac_, "bonmin.");
  options->GetNumericValue("maxmin_crit_no_sol", maxmin_crit_no_sol_, "bonmin.");
  options->GetNumericValue("maxmin_crit_have_sol", maxmin_crit_have_sol_, "bonmin.");

  /** Set values of standard branching options.*/
  int numberObjects = solver_->numberObjects();
  pseudoCosts_.initialize(numberObjects);
  int numberBeforeTrusted = b.getIntParameter(Bonmin::BabSetupBase::MinReliability);
  pseudoCosts_.setNumberBeforeTrusted(numberBeforeTrusted);

  setNumberStrong(b.getIntParameter(Bonmin::BabSetupBase::NumberStrong));

  /** Get values of options specific to CouenneChooseStrong.*/
  if (!options->GetIntegerValue("number_before_trust_list", numberBeforeTrustedList_, "bonmin.")) {
    // default is to use the same value as for numberBeforeTrusted
    numberBeforeTrustedList_ = numberBeforeTrusted;
  }
  options->GetIntegerValue("number_strong_branch_root", numberStrongRoot_, "bonmin.");

}

CouenneChooseStrong::CouenneChooseStrong(const CouenneChooseStrong & rhs) :
  OsiChooseStrong(rhs),
  cbc_model_(rhs.cbc_model_),
  only_pseudo_when_trusted_(rhs.only_pseudo_when_trusted_),
  maxmin_crit_no_sol_(rhs.maxmin_crit_no_sol_),
  maxmin_crit_have_sol_(rhs.maxmin_crit_have_sol_),
  setup_pseudo_frac_(rhs.setup_pseudo_frac_),
  numberBeforeTrustedList_(rhs.numberBeforeTrustedList_),
  numberStrongRoot_(rhs.numberStrongRoot_),
  problem_(rhs.problem_)
{
  jnlst_ = rhs.jnlst_;
  bb_log_level_ = rhs.bb_log_level_;
  DBG_ASSERT(IsValid(jnlst_));
}

CouenneChooseStrong &
CouenneChooseStrong::operator=(const CouenneChooseStrong & rhs)
{
  if (this != &rhs) {
    OsiChooseVariable::operator=(rhs);
    jnlst_ = rhs.jnlst_;
    bb_log_level_ = rhs.bb_log_level_;
    cbc_model_ = rhs.cbc_model_;
    only_pseudo_when_trusted_ = rhs.only_pseudo_when_trusted_;
    maxmin_crit_no_sol_ = rhs.maxmin_crit_no_sol_;
    maxmin_crit_have_sol_ = rhs.maxmin_crit_have_sol_;
    setup_pseudo_frac_ = rhs.setup_pseudo_frac_;
    numberBeforeTrustedList_ = rhs.numberBeforeTrustedList_;
    numberStrongRoot_ = rhs.numberStrongRoot_;
    problem_ = rhs.problem_;
  }
  return *this;
}

OsiChooseVariable *
CouenneChooseStrong::clone() const
{
  return new CouenneChooseStrong(*this);
}

CouenneChooseStrong::~CouenneChooseStrong ()
{
}

void
CouenneChooseStrong::computeMultipliers(double& upMult, double& downMult) const
{
  const double* upTotalChange = pseudoCosts_.upTotalChange();
  const double* downTotalChange = pseudoCosts_.downTotalChange();
  const int* upNumber = pseudoCosts_.upNumber();
  const int* downNumber = pseudoCosts_.downNumber();
  double sumUp=0.0;
  double numberUp=0.0;
  double sumDown=0.0;
  double numberDown=0.0;
  for (int i=numberObjects() - 1; i >= 0; --i) {
    sumUp += upTotalChange[i];
    numberUp += upNumber[i];
    sumDown += downTotalChange[i];
    numberDown += downNumber[i];
    if (bb_log_level_>4) {
      printf("%3d up %3d  %15.8e  down %3d  %15.8e\n",i, upNumber[i], upTotalChange[i], downNumber[i], downTotalChange[i]);
    }
  }
  upMult=(1.0+sumUp)/(1.0+numberUp);
  downMult=(1.0+sumDown)/(1.0+numberDown);

  if (bb_log_level_>4) {
    printf("upMultiplier = %e downMultiplier = %e\n", upMult, downMult);
  }
}

double
CouenneChooseStrong::computeUsefulness(const double MAXMIN_CRITERION,
				     const double upMult, const double downMult,
				     const double value,
				     const OsiObject* object, int i,
				     double& value2) const
{
  value2 = value;
  return value;
  double sumUp = pseudoCosts_.upTotalChange()[i]+1.0e-30;
  int numberUp = pseudoCosts_.upNumber()[i];
  double sumDown = pseudoCosts_.downTotalChange()[i]+1.0e-30;
  int numberDown = pseudoCosts_.downNumber()[i];
  double upEst = object->upEstimate();
  double downEst = object->downEstimate();
  upEst = numberUp ? ((upEst*sumUp)/numberUp) : (upEst*upMult);
  //if (numberUp<numberBeforeTrusted_)
  //  upEst *= (numberBeforeTrusted_+1.0)/(numberUp+1.0);
  downEst = numberDown ? ((downEst*sumDown)/numberDown) : (downEst*downMult);
  //if (numberDown<numberBeforeTrusted_)
  //  downEst *= (numberBeforeTrusted_+1.0)/(numberDown+1.0);
  double useful = ( MAXMIN_CRITERION*CoinMin(upEst,downEst) +
		    (1.0-MAXMIN_CRITERION)*CoinMax(upEst,downEst) );
  value2 = -COIN_DBL_MAX;
  if (numberUp   < numberBeforeTrustedList_ ||
      numberDown < numberBeforeTrustedList_) {
    value2 = value;
  }
  if (bb_log_level_>4) {
    printf("%3d value = %e upEstimate = %e downEstimate = %e infeas = %e value2 = %e\n", i,useful,upEst,downEst,value,value2);
  }
  return useful;
}

int
CouenneChooseStrong::setupList ( OsiBranchingInformation *info, bool initialize)
{
    problem_->update ((const CouNumber *) (info -> solution_), 
		      (const CouNumber *) (info -> lower_), 
		      (const CouNumber *) (info -> upper_));

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
  if (numberObjects>pseudoCosts_.numberObjects()) {
    //AW : How could that ever happen?  Right now, all old content is deleted!
    //    assert(false && "Right now, all old content is deleted!");
    // redo useful arrays
    pseudoCosts_.initialize(numberObjects);
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
  double upMultiplier, downMultiplier;
  computeMultipliers(upMultiplier, downMultiplier);

  // Say feasible
  bool feasible = true;
  const double MAXMIN_CRITERION = maxminCrit(info);
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
	double value2;
	value = computeUsefulness(MAXMIN_CRITERION,
				  upMultiplier, downMultiplier, value,
				  object[i], i, value2);
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
CouenneChooseStrong::chooseVariable(
  OsiSolverInterface * solver,
  OsiBranchingInformation *info,
  bool fixVariables)
{
  // We assume here that chooseVariable is called once at the very
  // beginning with fixVariables set to true.  This is then the root
  // node.
  bool isRoot = isRootNode(info);
  if (!isRoot) {
    numberStrongRoot_ = numberStrong_;
  }
  if (numberUnsatisfied_) {
    const double* upTotalChange = pseudoCosts_.upTotalChange();
    const double* downTotalChange = pseudoCosts_.downTotalChange();
    const int* upNumber = pseudoCosts_.upNumber();
    const int* downNumber = pseudoCosts_.downNumber();
    int numberBeforeTrusted = pseudoCosts_.numberBeforeTrusted();
    int numberLeft = CoinMin(numberStrongRoot_-numberStrongDone_,numberUnsatisfied_);
    int numberToDo=0;
    resetResults(numberLeft);
    int returnCode=0;
    bestObjectIndex_ = -1;
    bestWhichWay_ = -1;
    firstForcedObjectIndex_ = -1;
    firstForcedWhichWay_ =-1;
    double bestTrusted=-COIN_DBL_MAX;
    for (int i=0;i<numberLeft;i++) {
      int iObject = list_[i];
      if (numberBeforeTrusted<0||
	  (only_pseudo_when_trusted_ && number_not_trusted_>0) ||
	  !isRoot && (upNumber[iObject]<numberBeforeTrusted ||
		      downNumber[iObject]<numberBeforeTrusted )||
	  isRoot && (!upNumber[iObject] && !downNumber[iObject]) ) {
	results_[numberToDo++] = OsiHotInfo(solver, info,
					    solver->objects(), iObject);
      } else {
	const OsiObject * obj = solver->object(iObject);
	double upEstimate = (upTotalChange[iObject]*obj->upEstimate())/upNumber[iObject];
	double downEstimate = (downTotalChange[iObject]*obj->downEstimate())/downNumber[iObject];
	double MAXMIN_CRITERION = maxminCrit(info);
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
      returnCode = doStrongBranching(solver, info, numberToDo, 1);
      if (bb_log_level_>=3) {
	const char* stat_msg[] = {"NOTDON", "FEAS", "INFEAS", "NOFINI"};
	printf("BON0001I           DownStat    DownChange     UpStat      UpChange\n");
	for (int i = 0; i<numResults_; i++) {
	  double up_change = results_[i].upChange();
	  double down_change = results_[i].downChange();
	  int up_status = results_[i].upStatus();
	  int down_status = results_[i].downStatus();
	  printf("BON0002I    %3d    %6s    %13.6e   %6s    %13.6e\n",
		 i, stat_msg[down_status+1], down_change, stat_msg[up_status+1], up_change);
	}
      }
      if (returnCode>=0&&returnCode<=2) {
	if (returnCode) {
	  returnCode=4;
	  if (bestObjectIndex_>=0)
	    returnCode=3;
	}
	for (int i=0;i<numResults_;i++) {
	  int iObject = results_[i].whichObject();
	  double upEstimate;
	  if (results_[i].upStatus()!=1) {
	    assert (results_[i].upStatus()>=0);
	    upEstimate = results_[i].upChange();
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
	  if (results_[i].downStatus()!=1) {
	    assert (results_[i].downStatus()>=0);
	    downEstimate = results_[i].downChange();
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
	  double MAXMIN_CRITERION = maxminCrit(info);
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


bool CouenneChooseStrong::isRootNode(const OsiBranchingInformation *info) const {
  return info->depth_ == 0;
}

double
CouenneChooseStrong::maxminCrit(const OsiBranchingInformation *info) const {
  double retval = maxmin_crit_no_sol_;
#if 0
  if (cbc_model_) {
    // FIXME: should be replaced by info->stateOfSearch_
    const int stateOfSearch = cbc_model_->stateOfSearch();
    const int depth = info->depth_;
    if (stateOfSearch>1 && depth>10) {
      retval = maxmin_crit_have_sol_;
    }
  }
#endif
  return retval;
}

// Given a candidate  fill in useful information e.g. estimates
void 
CouenneChooseStrong::updateInformation(const OsiBranchingInformation *info,
				     int branch, OsiHotInfo * hotInfo)
{
  int index = hotInfo->whichObject();
  assert (index<solver_->numberObjects());
  const OsiObject * object = info->solver_->object(index);
  assert (object->upEstimate()>0.0&&object->downEstimate()>0.0);
  assert (branch<2);
  double* upTotalChange = pseudoCosts_.upTotalChange();
  double* downTotalChange = pseudoCosts_.downTotalChange();
  int* upNumber = pseudoCosts_.upNumber();
  int* downNumber = pseudoCosts_.downNumber();
   if (branch) {
    //if (hotInfo->upStatus()!=1) {
    // AW: Let's update the pseudo costs only if the strong branching
    // problem was marked as "solved"
    if (hotInfo->upStatus()==0) {
      assert (hotInfo->upStatus()>=0);
      upTotalChange[index] += hotInfo->upChange()/object->upEstimate();
      upNumber[index]++;
    } else if (hotInfo->upStatus()==1) {
      // infeasible - just say expensive
      upNumber[index]++;
      if (info->cutoff_<1.0e50)
	upTotalChange[index] += 2.0*(info->cutoff_-info->objectiveValue_)/object->upEstimate();
      else
	upTotalChange[index] += 2.0*fabs(info->objectiveValue_)/object->upEstimate();
    }
  } else {
    if (hotInfo->downStatus()==0) {
      assert (hotInfo->downStatus()>=0);
      downTotalChange[index] += hotInfo->downChange()/object->downEstimate();
      downNumber[index]++;
    } else if (hotInfo->upStatus()==1) {
      downNumber[index]++;
      // infeasible - just say expensive
      if (info->cutoff_<1.0e50)
	downTotalChange[index] += 2.0*(info->cutoff_-info->objectiveValue_)/object->downEstimate();
      else
	downTotalChange[index] += 2.0*fabs(info->objectiveValue_)/object->downEstimate();
    }
  }  
}

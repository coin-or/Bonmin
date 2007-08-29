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

namespace Bonmin {

BonChooseVariable::BonChooseVariable(OsiTMINLPInterface * solver) :
  OsiChooseStrong(solver),
  MAXMIN_CRITERION(0.85) 
{
  SmartPtr<TNLPSolver> tnlp_solver =
    static_cast<TNLPSolver *> (solver->solver());
  DBG_ASSERT(IsValid(tnlp_solver));
  jnlst_ = tnlp_solver->Jnlst();
  DBG_ASSERT(IsValid(jnlst_));
  SmartPtr<OptionsList> options = tnlp_solver->Options();

  options->GetIntegerValue("bb_log_level", bb_log_level_, "bonmin.");
}

BonChooseVariable::BonChooseVariable(const BonChooseVariable & rhs) :
  OsiChooseStrong(rhs),
  MAXMIN_CRITERION(rhs.MAXMIN_CRITERION) 
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
    MAXMIN_CRITERION = rhs.MAXMIN_CRITERION;
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
  if (initialize) {
    status_=-2;
    delete [] goodSolution_;
    bestObjectIndex_=-1;
    numberStrongDone_=0;
    numberStrongIterations_ = 0;
    numberStrongFixed_ = 0;
    goodSolution_ = NULL;
    goodObjectiveValue_ = COIN_DBL_MAX;
  }
  numberOnList_=0;
  numberUnsatisfied_=0;
  int numberObjects = solver_->numberObjects();
  if (numberObjects>numberObjects_) {
    // redo useful arrays
    delete [] upTotalChange_;
    delete [] downTotalChange_;
    delete [] upNumber_;
    delete [] downNumber_;
    numberObjects_ = solver_->numberObjects();
    upTotalChange_ = new double [numberObjects_];
    downTotalChange_ = new double [numberObjects_];
    upNumber_ = new int [numberObjects_];
    downNumber_ = new int [numberObjects_];
    CoinZeroN(upTotalChange_,numberObjects_);
    CoinZeroN(downTotalChange_,numberObjects_);
    CoinZeroN(upNumber_,numberObjects_);
    CoinZeroN(downNumber_,numberObjects_);
  }
  double check = -COIN_DBL_MAX;
  int checkIndex=0;
  int bestPriority=COIN_INT_MAX;
  int maximumStrong= CoinMin(numberStrong_,numberObjects) ;
  int putOther = numberObjects;
  int i;
  for (i=0;i<numberObjects;i++) {
    list_[i]=-1;
    useful_[i]=0.0;
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
  }
  double upMultiplier=(1.0+sumUp)/(1.0+numberUp);
  double downMultiplier=(1.0+sumDown)/(1.0+numberDown);
  // Say feasible
  bool feasible = true;
  for ( i=0;i<numberObjects;i++) {
    int way;
    double value = object[i]->infeasibility(info,way);
    if (value>0.0) {
      numberUnsatisfied_++;
      if (value==COIN_DBL_MAX) {
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
      } 
      if (priorityLevel==bestPriority) {
	// Modify value
	sumUp = upTotalChange_[i]+1.0e-30;
	numberUp = upNumber_[i];
	sumDown = downTotalChange_[i]+1.0e-30;
	numberDown = downNumber_[i];
	double upEstimate = object[i]->upEstimate();
	double downEstimate = object[i]->downEstimate();
	if (numberBeforeTrusted_>=0) {
	  if (shadowPossible<2) {
	    upEstimate = numberUp ? ((upEstimate*sumUp)/numberUp) : (upEstimate*upMultiplier);
	    if (numberUp<numberBeforeTrusted_)
	      upEstimate *= (numberBeforeTrusted_+1.0)/(numberUp+1.0);
	    downEstimate = numberDown ? ((downEstimate*sumDown)/numberDown) : (downEstimate*downMultiplier);
	    if (numberDown<numberBeforeTrusted_)
	      downEstimate *= (numberBeforeTrusted_+1.0)/(numberDown+1.0);
	  } else {
	    // use shadow prices always
	  }
	} else {
	  printf("numberBeforeTrusted_ < 0 is buggy\n");
	  abort();
	}
	value = MAXMIN_CRITERION*CoinMin(upEstimate,downEstimate) + (1.0-MAXMIN_CRITERION)*CoinMax(upEstimate,downEstimate);
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
      } else {
	// worse priority
	// to end
	assert (list_[putOther-1]<0);
	list_[--putOther]=i;
	maximumStrong = CoinMin(maximumStrong,putOther);
      }
    }
  }
  // Get list
  numberOnList_=0;
  if (feasible) {
    for (i=0;i<CoinMin(maximumStrong,putOther);i++) {
      if (list_[i]>=0) {
	list_[numberOnList_]=list_[i];
	useful_[numberOnList_++]=-useful_[i];
      }
    }
    if (numberOnList_) {
      // Sort 
      CoinSort_2(useful_,useful_+numberOnList_,list_);
      // move others
      i = numberOnList_;
      for (;putOther<numberObjects;putOther++) 
	list_[i++]=list_[putOther];
      assert (i==numberUnsatisfied_);
      if (!numberStrong_)
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
  if (bb_log_level_>4) {
    for (int i=0; i<Min(numberUnsatisfied_,numberStrong_); i++)
      printf("list_[%5d] = %5d, usefull_[%5d] = %23.16e\n", i,list_[i],i,useful_[i]);
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
  if (numberUnsatisfied_) {
    int numberLeft = CoinMin(numberStrong_-numberStrongDone_,numberUnsatisfied_);
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
      if (numberBeforeTrusted_<0||upNumber_[iObject]<numberBeforeTrusted_||downNumber_[iObject]<numberBeforeTrusted_) {
	results[numberToDo] = OsiHotInfo(solver,info,const_cast<const OsiObject **> (solver->objects()),iObject);
	temp[numberToDo++]=iObject;
      } else if (bestObjectIndex_<0) {
	const OsiObject * obj = solver->object(iObject);
	bestObjectIndex_=iObject;
	double upEstimate = (upTotalChange_[iObject]*obj->upEstimate())/upNumber_[iObject];
	double downEstimate = (downTotalChange_[iObject]*obj->downEstimate())/downNumber_[iObject];
	bestWhichWay_ = upEstimate>downEstimate ? 0 : 1;
	double value = MAXMIN_CRITERION*CoinMin(upEstimate,downEstimate) + (1.0-MAXMIN_CRITERION)*CoinMax(upEstimate,downEstimate);
	bestTrusted = value;
      }
    }
    int numberFixed=0;
    if (numberToDo) {
      int numberDone=0;
      returnCode = doBonStrongBranching(solver,info,numberToDo,results,1,numberDone);
      if (bb_log_level_>3) {
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
    } else {
#if 0
      // infeasible - just say expensive
      if (info->cutoff_<1.0e50)
	upTotalChange_[index] += 2.0*(info->cutoff_-info->objectiveValue_)/object->upEstimate();
      else
	upTotalChange_[index] += 2.0*fabs(info->objectiveValue_)/object->upEstimate();
#endif
    }
  } else {
    if (hotInfo->downStatus()==0) {
      assert (hotInfo->downStatus()>=0);
      downTotalChange_[index] += hotInfo->downChange()/object->downEstimate();
      downNumber_[index]++;
    } else {
#if 0
      // infeasible - just say expensive
      if (info->cutoff_<1.0e50)
	downTotalChange_[index] += 2.0*(info->cutoff_-info->objectiveValue_)/object->downEstimate();
      else
	downTotalChange_[index] += 2.0*fabs(info->objectiveValue_)/object->downEstimate();
#endif
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
  if (branch) {
    if (status!=1) {
      assert (status>=0);
      upTotalChange_[index] += changeInObjective/changeInValue;
      upNumber_[index]++;
    }
  } else {
    if (status!=1) {
      assert (status>=0);
      downTotalChange_[index] += changeInObjective/changeInValue;
      downNumber_[index]++;
    }
  }  
}

}/* Ends Bonmin's namespace.*/

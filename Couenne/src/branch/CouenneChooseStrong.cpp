/*
 * Name:    CouenneChooseStrong.cpp
 * Authors: Andreas Waechter, IBM Corp.
 * Purpose: Strong branching objects for Couenne
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneChooseStrong.hpp"
#include "CoinTime.hpp"

namespace Bonmin {

  CouenneChooseStrong::CouenneChooseStrong(BabSetupBase &b, CouenneProblem* p) :
    BonChooseVariable(b, b.continuousSolver()),
    problem_(p)
  {}


  CouenneChooseStrong::CouenneChooseStrong(const CouenneChooseStrong& rhs) :
    BonChooseVariable(rhs),
    problem_(rhs.problem_)
  {}

  CouenneChooseStrong::~CouenneChooseStrong()
  {}

  OsiChooseVariable *
  CouenneChooseStrong::clone() const
  {
    return new CouenneChooseStrong(*this);
  }

  CouenneChooseStrong&
  CouenneChooseStrong::operator=(const CouenneChooseStrong & rhs)
  {
    if (this != &rhs) {
      BonChooseVariable::operator=(rhs);
      problem_ = rhs.problem_;
    }
    return *this;
  }

  int
  CouenneChooseStrong::doStrongBranching( OsiSolverInterface * solver, 
					  OsiBranchingInformation *info,
					  int numberToDo, int returnCriterion)
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
    int returnCode=0;
    double timeStart = CoinCpuTime();
    int iDo = 0;
    for (;iDo<numberToDo;iDo++) {
      HotInfo * result = results_() + iDo;
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
        End of evaluation for this candidate object. Possibilities are:
        * Both sides below cutoff; this variable is a candidate for branching.
        * Both sides infeasible or above the objective cutoff: no further action
        here. Break from the evaluation loop and assume the node will be purged
        by the caller.
        * One side below cutoff: Install the branch (i.e., fix the variable). Possibly break
        from the evaluation loop and assume the node will be reoptimised by the
        caller.
      */
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
    if(iDo < numberToDo) iDo++;
    assert(iDo <= (int) results_.size());
    results_.resize(iDo);
    delete [] saveLower;
    delete [] saveUpper;
    // Delete the snapshot
    solver->unmarkHotStart();
    return returnCode;
  }

}

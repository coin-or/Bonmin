// Copyright (C) 2006, 2008 International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <climits>
#include "BonNWayChoose.hpp"
#include "BonNWayObject.hpp"
#include "CoinTime.hpp"

#ifndef NDEBUG
#define ASSERTED_CAST static_cast
#else
#define ASSERTED_CAST dynamic_cast
#endif
namespace Bonmin
{


  BonNWayChoose::BonNWayChoose(BabSetupBase &b, const OsiSolverInterface* solver):
      OsiChooseVariable(solver),
      br_depth_(0),
      bounds_(),
      unit_changes_()
  {

    Ipopt::SmartPtr<Ipopt::OptionsList> options = b.options();
    options->GetNumericValue("time_limit", time_limit_, b.prefix());
    options->GetIntegerValue("strong_branch_depth", br_depth_, b.prefix());

    /** Set values of standard branching options.*/
    int numberObjects = solver_->numberObjects();
    std::cout<<"Number objects "<<numberObjects<<std::endl;
    start_time_ = CoinCpuTime();
    OsiObject ** object = solver->objects();
    for (int i=0;i<numberObjects;i++) {
       BonNWayObject * nway = dynamic_cast<BonNWayObject *>(object[i]);
      if(!nway) continue;
      start_nway_ = i;
      break;
    }
    numberObjects -= start_nway_;
    std::cout<<"Number nway "<<numberObjects<<std::endl;

  }

  BonNWayChoose::BonNWayChoose(const BonNWayChoose & rhs) :
      OsiChooseVariable(rhs),
      br_depth_(rhs.br_depth_),
      time_limit_(rhs.time_limit_),
      start_time_(rhs.start_time_),
      start_nway_(rhs.start_nway_),
      bounds_(rhs.bounds_),
      unit_changes_(rhs.unit_changes_)
  {
  }

  BonNWayChoose &
  BonNWayChoose::operator=(const BonNWayChoose & rhs)
  {
    if (this != &rhs) {
      br_depth_ = rhs.br_depth_;
      time_limit_ = rhs.time_limit_;
      start_time_ = rhs.start_time_;
      start_nway_ = rhs.start_nway_;
      OsiChooseVariable::operator=(rhs);
      bounds_ = rhs.bounds_;
      unit_changes_ = rhs.unit_changes_;
    }
    return *this;
  }

  OsiChooseVariable *
  BonNWayChoose::clone() const
  {
    return new BonNWayChoose(*this);
  }

  BonNWayChoose::~BonNWayChoose ()
  {
  }

  void
  BonNWayChoose::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("NWay Strong branching setup", RegisteredOptions::BonminCategory);
    roptions->AddLowerBoundedIntegerOption("strong_branch_depth",
                                           "To which level do we perform strong-branching",
                                           0,0,
                                           "");

  }

  double
  BonNWayChoose::compute_usefulness(int objectIndex, const OsiBranchingInformation * info) const
  {
    int nwayIndex = objectIndex - start_nway_;
    const double * solution = info->solution_;
    double obj_val = info->objectiveValue_;
    const double * lower = info->lower_;
    const double * upper = info->upper_;
    double integerTolerance = info->integerTolerance_;
    double r_val = 1;

    BonNWayObject * nway = ASSERTED_CAST<BonNWayObject *>(info->solver_->objects()[objectIndex]);
    assert(nway);
    //BonNWayObject * nway = dynamic_cast<BonNWayObject *>(info->solver_->objects()[objectIndex]);
    size_t n = nway->numberMembers();
    const int * vars = nway->members();
    int n_pt = 0;

    return compute_usefulness(info, n, vars, bounds_[nwayIndex], unit_changes_[nwayIndex]);
    for(size_t  i = 0 ; i < n ; i++){
      int iCol = vars[i];
      if(fabs(lower[iCol] - upper[iCol]) < integerTolerance) continue; //Variable is fixed
      assert(lower[iCol] < upper[iCol]);
      double residual = upper[iCol] - solution[iCol];
      r_val*=std::max(residual*unit_changes_[nwayIndex][i],bounds_[nwayIndex][i] - obj_val);
      n_pt++;
    } 
    //printf("r_val %g n_pt %i\n", r_val, n_pt);
    return r_val;
    return pow(r_val, 1./(double) n_pt);
  }

  double
  BonNWayChoose::compute_usefulness(const OsiBranchingInformation * info,
                size_t n, const int * vars,const std::vector<double> &bounds, const std::vector<double> &unit_changes) const
  {
    const double * solution = info->solution_;
    double obj_val = info->objectiveValue_;
    const double * lower = info->lower_;
    const double * upper = info->upper_;
    double integerTolerance = info->integerTolerance_;
    double r_val = 1;

    int n_pt = 0;

    for(size_t  i = 0 ; i < n ; i++){
      int iCol = vars[i];
      if(fabs(lower[iCol] - upper[iCol]) < integerTolerance) continue; //Variable is fixed
      assert(lower[iCol] < upper[iCol]);
      double residual = upper[iCol] - solution[iCol];
      r_val*=std::max(residual*unit_changes[i],bounds[i] - obj_val);
      n_pt++;
    } 
    //printf("r_val %g n_pt %i\n", r_val, n_pt);
    return r_val;
    return pow(r_val, 1./(double) n_pt)/n_pt;
  }
  int
  BonNWayChoose::setupList ( OsiBranchingInformation *info, bool initialize)
  {
    assert(initialize);
    if (initialize) {
      status_=-2;
      delete [] goodSolution_;
      bestObjectIndex_=-1;
      goodSolution_ = NULL;
      goodObjectiveValue_ = COIN_DBL_MAX;
    }
#if 0
    else {
      std::cerr<<"BonNWayChoose::setupList , Should not be called"
                <<" with initialize==false"<<std::endl;
      throw CoinError("BonNWayChoose","setupList","Should not be called with initialize==false");
    }
#endif


    numberOnList_=0;
    numberUnsatisfied_=0;
    int numberObjects = info->solver_->numberObjects() - start_nway_;
    int maximumStrong = 0;
    if(info->depth_ == 0){
      bounds_.resize(0);
      unit_changes_.resize(0);
      bounds_.resize(numberObjects, std::vector<double>(numberObjects,0.));
      unit_changes_.resize(numberObjects, std::vector<double>(numberObjects,0.));
    }
    else {
      assert(unit_changes_.size() == numberObjects);
      assert(bounds_.size() == unit_changes_.size());
    }
    if(info->depth_ > br_depth_)
      maximumStrong = 1;
    else
      maximumStrong = numberObjects;



    double check = -COIN_DBL_MAX;
    int checkIndex=0;
    int bestPriority=COIN_INT_MAX;
    int putOther = numberObjects;
    int i;
    for (i=0;i<numberObjects;i++) {
      list_[i]=-1;
      useful_[i]=0.0;
    }

    OsiObject ** object = info->solver_->objects();
    object += start_nway_;

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
          bestPriority = priorityLevel;
          check=-COIN_DBL_MAX;
          checkIndex=0;
        }
        if (priorityLevel==bestPriority) {
          // Modify value
          //Compute usefullness
          if(info->depth_ != 0){
            value = compute_usefulness(i + start_nway_, info); 
            //printf("Usefullness on pseudo cost %i : %g\n", i, value);
          }

          if (value>check) {
            //add to list
            int iObject = list_[checkIndex];
            if (iObject>=0) {
              assert (list_[putOther-1]<0);
              list_[--putOther]=iObject;  // to end
            }
            list_[checkIndex]= i + start_nway_;
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
               }
               else {
                 check=0.0;
                 checkIndex = j;
                 break;
               }
             }
          }
          else {
            // to end
            assert (list_[putOther-1]<0);
            list_[--putOther]=i + start_nway_;
            maximumStrong = CoinMin(maximumStrong,putOther);
          }
        }
        else {
          // worse priority
          // to end
          assert (list_[putOther-1]<0);
          list_[--putOther]=i + start_nway_;
          maximumStrong = CoinMin(maximumStrong,putOther);
        }
      }
    }

    // Get list
    numberOnList_=0;
    if (feasible) {
      maximumStrong = CoinMin(maximumStrong,putOther);
      for (i=0;i<maximumStrong;i++) {
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
    }
    else {
      // not feasible
      numberUnsatisfied_=-1;
    }
    // Get rid of any shadow prices info
    info->defaultDual_ = -1.0; // switch off
    delete [] info->usefulRegion_;
    delete [] info->indexRegion_;

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
  BonNWayChoose::chooseVariable(
    OsiSolverInterface * solver,
    OsiBranchingInformation *info,
    bool fixVariables)
  {
    if(!numberUnsatisfied_) return 1;//Node is feasible



    if(info->depth_ > br_depth_){
      const double * solution = info->solution_;
      double obj_val = info->objectiveValue_;
      double cutoff = info->cutoff_;
      const double * lower = info->lower_;
      const double * upper = info->upper_;
      double integerTolerance = info->integerTolerance_;


      // See if quick variable fixing can be done
      int n_fixed = 0;
      for(int i = 0 ; i < numberUnsatisfied_ ; i++){
         int iObject = list_[i];
         int nwayIndex = iObject - start_nway_;
         const BonNWayObject * nway = ASSERTED_CAST<const BonNWayObject *>(solver->object(iObject));

         size_t n = nway->numberMembers();
         const int * vars = nway->members();
         for(size_t  j = 0 ; j < n ; j++){
           int iCol = vars[j];
           if(upper[iCol] < lower[iCol] + 0.5) continue;//variable already fixed to lower buund
           if(bounds_[nwayIndex][j] > cutoff){//It can be fixed
              solver->setColUpper(iCol, lower[iCol]);
              n_fixed ++;
           }
         }
      }
      if(n_fixed)
        printf("Fixed %i variables\n", n_fixed);

      assert(bounds_.size() == unit_changes_.size());
      assert(unit_changes_.size() == info->solver_->numberObjects() - start_nway_);
      //Just take the most usefull
      bestObjectIndex_ = list_[0];
      bestWhichWay_ = 1; 
      OsiObject * obj = solver->objects()[bestObjectIndex_];
      obj->setWhichWay(bestWhichWay_);

#define OUT_BRANCH
#ifdef OUT_BRANCH
      int nwayIndex = bestObjectIndex_ - start_nway_;
      BonNWayObject * nway = ASSERTED_CAST<BonNWayObject *>
                        (info->solver_->objects()[bestObjectIndex_]);
      assert(nway);
      size_t n = nway->numberMembers();
      const int * vars = nway->members();
      int n_pt = 0;
      double r_val = 1; 
      for(size_t  i = 0 ; i < n ; i++){
        int iCol = vars[i];
        if(fabs(lower[iCol] - upper[iCol]) < integerTolerance) continue; //Variable is fixed
        assert(lower[iCol] < upper[iCol]);
        double residual = upper[iCol] - solution[iCol];
        r_val*=std::max(residual*unit_changes_[nwayIndex][i],bounds_[nwayIndex][i] - obj_val);
        if(bounds_[nwayIndex][i] < cutoff) n_pt++;
      }
      printf("level %i branch on %i: %i branches, bound %g usefullness %g.\n",
             info->depth_, bestObjectIndex_, n_pt, obj_val, r_val);
#endif 

      if(n_fixed) return 2;
      return 0;
     }


    printf("Restarting strong branching loop....\n\n");

    int numberLeft = numberUnsatisfied_;//Always do full strong at root
    int returnCode=0;
    bestObjectIndex_ = -1;
    bestWhichWay_ = -1;
    firstForcedObjectIndex_ = -1;
    firstForcedWhichWay_ =-1;
    double best_score = -COIN_DBL_MAX;
    int bestPriority=0;

    int n = solver->getNumCols();
    std::vector<double> saveLower(n);
    std::vector<double> saveUpper(n);
    std::copy(info->lower_, info->lower_ + n, saveLower.begin());
    std::copy(info->upper_, info->upper_ + n, saveUpper.begin());


    solver->markHotStart();
    for (int i=0;i<numberLeft;i++) {
      int iObject = list_[i];
      const int objectPriority = solver->object(iObject)->priority();
      if (objectPriority >= bestPriority){
             bestPriority = objectPriority;
      }
      else break;
      double score;
      int r_val = doStrongBranching(solver, info, iObject, saveLower.data(),
                                     saveUpper.data(), score);
      if(r_val == -1) {
         std::cerr<< "This is Infeasible"<<std::endl;
         returnCode = -1;
         break;
      }

       if (r_val == 1 && fixVariables && 0) {//for now don't as this is very unlikely
             //and would require nontrivial coding.
             const OsiObject * obj = solver->object(iObject);
             OsiBranchingObject * branch = obj->createBranch(solver,info,0);
             branch->branch(solver);
             delete branch;
      }
      if(r_val == 1 && info->depth_ == 0) returnCode=2;
      //Compute Score
      printf("Usefullness from strong branching on %i : %g\n", iObject - start_nway_, score);
      if(score > best_score){//We have a winner
        best_score = score;
        bestObjectIndex_ = iObject;
        bestWhichWay_ = 0;
      }
      if (r_val==3) {
          returnCode = 3;
          // max time - just choose one
          if(bestObjectIndex_ < 0){
            bestObjectIndex_ = list_[0];
            bestWhichWay_ = 0;
          }
          break;
        }
    }
    solver->unmarkHotStart();
    assert(returnCode == 0);
    return returnCode;
  }

  /*  This is a utility function which does strong branching on
      one nway object and stores the results in appropriate arrays of the class
      and maybe more.
      It returns -
      -1 - branch was infeasible both ways
       0 - nothing can be fixed
       1 - object can be fixed (returnCriterion==0)
       3 - time limit
  */
  int 
  BonNWayChoose::doStrongBranching( OsiSolverInterface * solver, 
  				    OsiBranchingInformation *info,
  				    int objectIndex,
                                    double * saveLower,
                                    double * saveUpper, double & score)
  {
    int nwayIndex = objectIndex - start_nway_;
    const double * lower = info->lower_;
    const double * upper = info->upper_;

    int numberColumns = solver->getNumCols();
    double timeStart = CoinCpuTime();

    int numberObjects = info->solver_->numberObjects();
    const BonNWayObject * nway = ASSERTED_CAST<const BonNWayObject *>(solver->object(objectIndex));
    assert(nway);
    BonNWayBranchingObject * branch = ASSERTED_CAST<BonNWayBranchingObject *>(nway->createBranch(solver, info, 1));

    int branches_left = branch->numberBranchesLeft();
    int number_branches = branch->numberBranchesLeft();
    int n_can_be_fixed = 0;

    double big_val = 3*info->cutoff_;// For infeasibles
    if(big_val > 1e10){ big_val = 10*info->objectiveValue_;}
    big_val += fabs(big_val)*1e-5;
    std::vector<double> unit_changes(numberObjects - start_nway_, 0.);
    while(branches_left){ 

      branch->branch(solver);
      int v_br = branch->var_branched_on();
      int s_br = branch->seq_branched_on();
      double residual = upper[v_br] - info->solution_[v_br];
      solver->solveFromHotStart() ;
      numberStrongIterations_ += solver->getIterationCount();
      numberStrongDone_++;

      double obj_val = solver->getObjValue();

      if(solver->isProvenPrimalInfeasible() ||
         (solver->isProvenOptimal() && obj_val > info->cutoff_)){//infeasible
         n_can_be_fixed++;
         if(info->depth_ == 0){
         bounds_[nwayIndex][s_br] = big_val;
         if(0 && obj_val > 1e10){//failure be less violent
           bounds_[nwayIndex][s_br] = - 1e50;
         }
         //Fix variable
         }
         unit_changes[s_br] = (big_val - info->objectiveValue_)/residual;
         printf("Fixing variable %i to 0\n", v_br);
         saveUpper[v_br] = saveLower[v_br];
      }
      else{
         if(info->depth_ == 0){
          bounds_[nwayIndex][s_br] = obj_val;
         }
         unit_changes[s_br] = (obj_val - info->objectiveValue_)/residual;
      }

      // Restore bounds
      for (int j=0;j<numberColumns;j++) {
        if (saveLower[j] != lower[j])
  	solver->setColLower(j,saveLower[j]);
        if (saveUpper[j] != upper[j])
  	solver->setColUpper(j,saveUpper[j]);
      }
      branches_left = branch->numberBranchesLeft();
    }

    score = compute_usefulness(info, nway->numberMembers(), nway->members(), bounds_[nwayIndex], unit_changes);
    if(info->depth_ > 0) score = pow(score, 1./(number_branches-n_can_be_fixed));
    if(info->depth_ == 0){//At root bounds contains valid bound on obj after branching, remember
        nway->set_bounds(bounds_[nwayIndex]);
        unit_changes_[nwayIndex] = unit_changes;
    }

    if(n_can_be_fixed == number_branches){
       return -1;
    }
    //if(n_can_be_fixed == number_branches - 1){
    if(n_can_be_fixed){
       return 1;
    }
    bool hitMaxTime = ( CoinCpuTime()-timeStart > info->timeRemaining_)
                        || ( CoinCpuTime() - start_time_ > time_limit_);
   if (hitMaxTime) {
        return 3;
    }
   return 0;
  }

}/* Ends Bonmin's namespace.*/


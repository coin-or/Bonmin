// Copyright (C) 2006, 2008 International Business Machines
// Corporation and others.  All Rights Reserved.

#include <climits>
#include "CoinPragma.hpp"
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
      unit_changes_(),
      num_ps_costs_(),
      num_eval_(),
      geo_means_(0)
  {
    Ipopt::SmartPtr<Ipopt::OptionsList> options = b.options();
    options->GetNumericValue("time_limit", time_limit_, b.prefix());
    options->GetNumericValue("cutoff_multiplier", cutoff_multiplier_, b.prefix());
    options->GetNumericValue("pseudocost_trust_value", pseudocost_trust_value_, b.prefix());
    options->GetIntegerValue("strong_branch_depth", br_depth_, b.prefix());
    options->GetIntegerValue("nway_branch_log_level", log_, b.prefix());
    options->GetEnumValue("do_fixings", do_fixings_, b.prefix());
    options->GetEnumValue("use_geo_means", geo_means_, b.prefix());
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
  }

  BonNWayChoose::BonNWayChoose(const BonNWayChoose & rhs) :
      OsiChooseVariable(rhs),
      br_depth_(rhs.br_depth_),
      do_fixings_(rhs.do_fixings_),
      cutoff_multiplier_(rhs.cutoff_multiplier_),
      pseudocost_trust_value_(rhs.pseudocost_trust_value_),
      time_limit_(rhs.time_limit_),
      start_time_(rhs.start_time_),
      start_nway_(rhs.start_nway_),
      log_(rhs.log_),
      bounds_(rhs.bounds_),
      unit_changes_(rhs.unit_changes_),
      num_ps_costs_(rhs.num_ps_costs_),
      num_eval_(rhs.num_eval_),
      geo_means_(rhs.geo_means_)
  {
  }

  BonNWayChoose &
  BonNWayChoose::operator=(const BonNWayChoose & rhs)
  {
    if (this != &rhs) {
      br_depth_ = rhs.br_depth_;
      do_fixings_ = rhs.do_fixings_;
      cutoff_multiplier_ = rhs.cutoff_multiplier_;
      pseudocost_trust_value_ = rhs.pseudocost_trust_value_;
      time_limit_ = rhs.time_limit_;
      start_time_ = rhs.start_time_;
      log_ = rhs.log_;
      start_nway_ = rhs.start_nway_;
      OsiChooseVariable::operator=(rhs);
      bounds_ = rhs.bounds_;
      unit_changes_ = rhs.unit_changes_;
      num_ps_costs_ = rhs.num_ps_costs_;
      num_eval_ = rhs.num_eval_;
      geo_means_ = rhs.geo_means_;
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
    roptions->AddLowerBoundedIntegerOption("nway_branch_log_level",
                                           "Log level for the branching on nways",
                                           0,1,
                                           "");

    roptions->AddLowerBoundedIntegerOption("strong_branch_depth",
                                           "To which level do we perform strong-branching",
                                           0,0,
                                           "");

    roptions->AddLowerBoundedNumberOption("cutoff_multiplier",
                                           "multiplier applied to cutoff_ for computing pseudo-cost of infeasible sub-problems",
                                           1.,0,3.,
                                           "");

    roptions->AddLowerBoundedNumberOption("pseudocost_trust_value",
                                           "Trust pseudo cost of best nway object if it is above this value",
                                           0.,0,0,
                                           "");

    roptions->AddStringOption2("use_geo_means", "Use geometrical means to average pseudo-costs",
                               "yes", 
                               "no", "Use artihmetical means",
                               "yes", "Use geometrical means","");

    roptions->AddStringOption4("do_fixings",
        "Do we fix variables in strong branching?",
        "all",
        "none", "Don't do any.",
        "in-tree", "Fix only variables in the tree",
        "strong-branching", "Fix variable in strong branching only",
        "all", "Fix whenever possible",
        "");


  }

  double
  BonNWayChoose::compute_usefulness(int objectIndex, const OsiBranchingInformation * info) const
  {
    int nwayIndex = objectIndex - start_nway_;

    BonNWayObject * nway = ASSERTED_CAST<BonNWayObject *>(info->solver_->objects()[objectIndex]);
    assert(nway);
    size_t n = nway->numberMembers();
    const int * vars = nway->members();

    std::vector<double> unit_changes(unit_changes_[nwayIndex]);
    if(geo_means_){
      for(size_t k = 0 ; k < unit_changes.size() ; k++) unit_changes[k] = (num_ps_costs_[nwayIndex][k]) ? pow(unit_changes[k], 1./(double) num_ps_costs_[nwayIndex][k]): unit_changes[k];
    }
    else {
    for(size_t k = 0 ; k < unit_changes.size() ; k++) unit_changes[k] = (num_ps_costs_[nwayIndex][k]) ? unit_changes[k]/(double) num_ps_costs_[nwayIndex][k]: unit_changes[k];
    }


    double r_val = compute_usefulness(info, n, vars, bounds_[nwayIndex], unit_changes);
    return r_val;
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
    double cutoff = info->cutoff_*cutoff_multiplier_;
    double r_val = (geo_means_) ? 1 : 0;

    for(size_t  i = 0 ; i < n ; i++){
      int iCol = vars[i];
      if(fabs(lower[iCol] - upper[iCol]) < integerTolerance) {
        //r_val *= (cutoff - obj_val);
        continue; //Variable is fixed
      }
      assert(lower[iCol] < upper[iCol]);
      double residual = upper[iCol] - solution[iCol];
      double score = std::min(cutoff - obj_val, std::max(residual*unit_changes[i],bounds[i] - obj_val));
      if(geo_means_)
        r_val*=score;
      else
        r_val += score;
    } 
    return r_val;
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


    numberOnList_=0;
    numberUnsatisfied_=0;
    int numberObjects = info->solver_->numberObjects() - start_nway_;
    double cutoff = info->cutoff_;
    if(info->depth_ == 0){
      bounds_.resize(0);
      unit_changes_.resize(0);
      bounds_.resize(numberObjects, std::vector<double>(numberObjects,cutoff));
      unit_changes_.resize(numberObjects, std::vector<double>(numberObjects,0));
      num_eval_.resize(numberObjects, 0);
      num_ps_costs_.resize(numberObjects, std::vector<int>(numberObjects,0));
    }
    else {
      assert(unit_changes_.size() == numberObjects);
      assert(bounds_.size() == unit_changes_.size());
    }
    //Always allow all objects on the list  
    int maximumStrong = numberObjects;

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


    double cutoff = info->cutoff_;
    double obj_val = info->objectiveValue_;
    if(info->depth_ > br_depth_ &&  ( (- useful_[0]   > pseudocost_trust_value_ )|| num_eval_[list_[0] - start_nway_] >= 18) ){
      const double * lower = info->lower_;
      const double * upper = info->upper_;


      // See if quick variable fixing can be done
      int n_fixed = 0;
      if(do_fixings_ > 1){
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
        if(n_fixed && log_ > 1)
          printf("NWAY: Fixed %i variables\n", n_fixed);
      }

      assert(bounds_.size() == unit_changes_.size());
      assert(unit_changes_.size() == info->solver_->numberObjects() - start_nway_);
      //Just take the most usefull
      bestObjectIndex_ = list_[0];
      bestWhichWay_ = 1; 
      OsiObject * obj = solver->objects()[bestObjectIndex_];
      obj->setWhichWay(bestWhichWay_);

      if(log_ > 1){
        printf("level %i branch on %i bound %g usefullness %g.\n",
               info->depth_, bestObjectIndex_ - start_nway_, obj_val, - useful_[0]);
      }
      if(n_fixed) return 2;
      return 0;
     }


    if(log_ > 0)
      printf("Restarting strong branching loop....\n\n");

    numberStrongIterations_ = 0;
    numberStrongDone_ = 0;
    int numberLeft = numberOnList_;//Always do full strong at root
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
         if(log_ > 0)
           std::cout<<"This is Infeasible"<<std::endl;
         returnCode = -1;
         break;
      }

      if(do_fixings_ > 1 && r_val == 1 && info->depth_ == 0) returnCode=2;
      //Compute Score
      if(log_ > 0)
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

    double big_val = cutoff_multiplier_*info->cutoff_;// For infeasibles
    if(big_val > 1e10){ big_val = 10*info->objectiveValue_;}
    big_val += fabs(big_val)*1e-5;
    std::vector<double> unit_changes(numberObjects - start_nway_, -DBL_MAX);
    //std::vector<double> unit_changes(numberObjects - start_nway_, 0);
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
         if(info->depth_ == 0){
         bounds_[nwayIndex][s_br] = big_val;
         }
         unit_changes[s_br] = (big_val - info->objectiveValue_)/residual;
         if(do_fixings_ > 1){
           n_can_be_fixed++;
           if(log_ > 0)
           printf("Fixing variable %i to 0 the cutoff is %g\n", v_br, big_val);
           saveUpper[v_br] = saveLower[v_br];
         }
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
    if(info->depth_ == 0){//At root bounds contains valid bound on obj after branching, remember
        if(do_fixings_ == 1 || do_fixings_ == 3)
        nway->set_bounds(bounds_[nwayIndex]);
        for(size_t k = 0 ; k < unit_changes.size() ; k++){
           num_ps_costs_[nwayIndex][k]=1;
        }
        unit_changes_[nwayIndex] = unit_changes;
        num_eval_[nwayIndex] = 1;
    }
    else if (n_can_be_fixed < number_branches -1){
       num_eval_[nwayIndex]++;
       for(size_t k = 0 ; k < unit_changes.size() ; k++){
         if(unit_changes[k] > 0.){
           if(geo_means_)
              unit_changes_[nwayIndex][k] *= unit_changes[k];
           else
              unit_changes_[nwayIndex][k] += unit_changes[k];
           num_ps_costs_[nwayIndex][k]++;
         }
       }
    }
    if(n_can_be_fixed == number_branches){
       return -1;
    }
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


// Copyright (C) 2006, 2008 International Business Machines
// Corporation and others.  All Rights Reserved.

#include <climits>
#include "CoinPragma.hpp"
#include "BonChooseVariable.hpp"
#include "CoinTime.hpp"
#include "IpBlas.hpp"
#include "BonMsgUtils.hpp"

// This couples Cbc code into Bonmin code...
#include "CbcModel.hpp"

namespace Bonmin
{

  BonChooseVariable::Messages::Messages():
      CoinMessages((int) BON_CHOOSE_MESSAGES_DUMMY_END)
  {
    strcpy(source_,"BON");
    ADD_MSG(PS_COST_HISTORY,std_m,6,"%3d up %3d  %.8e  down %3d  %.8e");
    ADD_MSG(PS_COST_MULT,std_m, 6, "upMultiplier = %e downMultiplier = %e");
    ADD_MSG(PS_COST_ESTIMATES, std_m, 6, "%3d value = %e upEstimate = %e downEstimate = %e infeas = %e value2 = %e");
    ADD_MSG(CANDIDATE_LIST,std_m,6,
        "list_[%5d] = %5d, usefull_[%5d] = %.16e %.16e");
    ADD_MSG(CANDIDATE_LIST2, std_m, 6,
        "list_[%3d] = %3d useful_[%3d] = %e");
    ADD_MSG(CANDIDATE_LIST3, std_m, 6,
        "list2[%3d] = %3d useful2[%3d] = %e");
    ADD_MSG(SB_START, std_m,5,
        " Starting strong branching. Obj. val = %g\n");
    ADD_MSG(SB_HEADER, std_m,5,
        "           Var    Value            DownStat    DownChange     UpStat      UpChange");
    ADD_MSG(SB_RES, std_m, 5,
        "    %3d    %3d    %.6e      %6s    %.6e   %6s    %.6e");
    ADD_MSG(BRANCH_VAR, std_m, 4, "Branched on variable %i, bestWhichWay: %i");
    ADD_MSG(CHOSEN_VAR, std_m, 4,"           Choosing %d");
    ADD_MSG(UPDATE_PS_COST, std_m, 4,"update %3d %3d %e %e %3d");
  }
  const std::string BonChooseVariable::CNAME = "BonChooseVariable";

  BonChooseVariable::BonChooseVariable(BabSetupBase &b, const OsiSolverInterface* solver):
      OsiChooseVariable(solver),
      results_(),
      cbc_model_(NULL),
      only_pseudo_when_trusted_(false),
      pseudoCosts_()
  {
    jnlst_ = b.journalist();
    Ipopt::SmartPtr<Ipopt::OptionsList> options = b.options();

    handler_ = new CoinMessageHandler;

    options->GetIntegerValue("bb_log_level", bb_log_level_, b.prefix());
    handler_->setLogLevel(bb_log_level_);
    options->GetNumericValue("time_limit", time_limit_, b.prefix());
    options->GetNumericValue("setup_pseudo_frac", setup_pseudo_frac_, b.prefix());
    options->GetNumericValue("maxmin_crit_no_sol", maxmin_crit_no_sol_, b.prefix());
    options->GetNumericValue("maxmin_crit_have_sol", maxmin_crit_have_sol_, b.prefix());
    options->GetEnumValue("trust_strong_branching_for_pseudo_cost",trustStrongForPseudoCosts_ , b.prefix());
    int sortCrit;
    options->GetEnumValue("candidate_sort_criterion", sortCrit, b.prefix());
#ifndef OLD_USEFULLNESS
    sortCrit_ = (CandidateSortCriterion) sortCrit;
#endif
    /** Set values of standard branching options.*/
    int numberObjects = solver_->numberObjects();
    //std::cout<<"Number objects "<<numberObjects<<std::endl;
    pseudoCosts_.initialize(numberObjects);
    int numberBeforeTrusted = b.getIntParameter(BabSetupBase::MinReliability);
    pseudoCosts_.setNumberBeforeTrusted(numberBeforeTrusted);

    setNumberStrong(b.getIntParameter(BabSetupBase::NumberStrong));

    /** Get values of options specific to BonChooseVariable.*/
    if (!options->GetIntegerValue("number_before_trust_list", numberBeforeTrustedList_, b.prefix())) {
      // default is to use the same value as for numberBeforeTrusted
      numberBeforeTrustedList_ = numberBeforeTrusted;
    }
    options->GetIntegerValue("number_strong_branch_root", numberStrongRoot_, b.prefix());
    options->GetIntegerValue("min_number_strong_branch", minNumberStrongBranch_, b.prefix());
    options->GetIntegerValue("number_look_ahead", numberLookAhead_, b.prefix());

    start_time_ = CoinCpuTime();
  }

  BonChooseVariable::BonChooseVariable(const BonChooseVariable & rhs) :
      OsiChooseVariable(rhs),
      results_(rhs.results_),
      time_limit_(rhs.time_limit_),
      start_time_(CoinCpuTime()),
      cbc_model_(rhs.cbc_model_),
      only_pseudo_when_trusted_(rhs.only_pseudo_when_trusted_),
      maxmin_crit_no_sol_(rhs.maxmin_crit_no_sol_),
      maxmin_crit_have_sol_(rhs.maxmin_crit_have_sol_),
      setup_pseudo_frac_(rhs.setup_pseudo_frac_),
      numberBeforeTrustedList_(rhs.numberBeforeTrustedList_),
      numberStrongRoot_(rhs.numberStrongRoot_),
#ifndef OLD_USEFULLNESS
      sortCrit_(rhs.sortCrit_),
#endif
      numberLookAhead_(rhs.numberLookAhead_),
      minNumberStrongBranch_(rhs.minNumberStrongBranch_),
      pseudoCosts_(rhs.pseudoCosts_),
      trustStrongForPseudoCosts_(rhs.trustStrongForPseudoCosts_)
  {
    jnlst_ = rhs.jnlst_;
    handler_ = rhs.handler_->clone();
    bb_log_level_ = rhs.bb_log_level_;
    DBG_ASSERT(IsValid(jnlst_));
  }

  BonChooseVariable &
  BonChooseVariable::operator=(const BonChooseVariable & rhs)
  {
    if (this != &rhs) {
      OsiChooseVariable::operator=(rhs);
      delete handler_;
      handler_ = rhs.handler_->clone();
      jnlst_ = rhs.jnlst_;
      bb_log_level_ = rhs.bb_log_level_;
      cbc_model_ = rhs.cbc_model_;
      only_pseudo_when_trusted_ = rhs.only_pseudo_when_trusted_;
      maxmin_crit_no_sol_ = rhs.maxmin_crit_no_sol_;
      maxmin_crit_have_sol_ = rhs.maxmin_crit_have_sol_;
      setup_pseudo_frac_ = rhs.setup_pseudo_frac_;
      numberBeforeTrustedList_ = rhs.numberBeforeTrustedList_;
      numberStrongRoot_ = rhs.numberStrongRoot_;
#ifndef OLD_USEFULLNESS
      sortCrit_ = rhs.sortCrit_;
#endif
      minNumberStrongBranch_ = rhs.minNumberStrongBranch_;
      pseudoCosts_ = rhs.pseudoCosts_;
      trustStrongForPseudoCosts_ = rhs.trustStrongForPseudoCosts_;
      numberLookAhead_ = rhs.numberLookAhead_;
      results_ = rhs.results_;
    }
    return *this;
  }

  OsiChooseVariable *
  BonChooseVariable::clone() const
  {
    return new BonChooseVariable(*this);
  }

  BonChooseVariable::~BonChooseVariable ()
  {
    delete handler_;
  }

  void
  BonChooseVariable::registerOptions(
    Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("Strong branching setup", RegisteredOptions::BonminCategory);
    roptions->AddStringOption4("candidate_sort_criterion",
        "Choice of the criterion to choose candidates in strong-branching",
        "best-ps-cost",
        "best-ps-cost","Sort by decreasing pseudo-cost",
        "worst-ps-cost", "Sort by increasing pseudo-cost",
        "most-fractional", "Sort by decreasing integer infeasibility",
        "least-fractional", "Sort by increasing integer infeasibility","");

    roptions->setOptionExtraInfo("candidate_sort_criterion",63);

    roptions->AddBoundedNumberOption("setup_pseudo_frac", "Proportion of strong branching list that has to be taken from most-integer-infeasible list.",
        0., false, 1., false, 0.5);
    roptions->setOptionExtraInfo("setup_pseudo_frac",63);
    roptions->AddBoundedNumberOption("maxmin_crit_no_sol", "Weight towards minimum in of lower and upper branching estimates when no solution has been found yet.",
        0., false, 1., false, 0.7);
    roptions->setOptionExtraInfo("maxmin_crit_no_sol",63);
    roptions->AddBoundedNumberOption("maxmin_crit_have_sol", "Weight towards minimum in of lower and upper branching estimates when a solution has been found.",
        0., false, 1., false, 0.1);
    roptions->setOptionExtraInfo("maxmin_crit_have_sol",63);
    roptions->AddLowerBoundedIntegerOption("number_before_trust_list",
        "Set the number of branches on a variable before its pseudo costs are to be believed during setup of strong branching candidate list.",
        -1, 0, "The default value is that of \"number_before_trust\"");
    roptions->setOptionExtraInfo("number_before_trust_list",63);
    roptions->AddLowerBoundedIntegerOption("number_strong_branch_root",
        "Maximum number of variables considered for strong branching in root node.",
        0, COIN_INT_MAX, "");
    roptions->setOptionExtraInfo("number_strong_branch_root",63);

    roptions->AddLowerBoundedIntegerOption("min_number_strong_branch", "Sets minimum number of variables for strong branching (overriding trust)",
        0, 0,"");
    roptions->setOptionExtraInfo("min_number_strong_branch",63);
    roptions->AddStringOption2("trust_strong_branching_for_pseudo_cost",
                               "Whether or not to trust strong branching results for updating pseudo costs.",
                               "yes",
                               "no","",
                               "yes","",
                               ""
                               );
    roptions->setOptionExtraInfo("trust_strong_branching_for_pseudo_cost", 63);

    roptions->AddLowerBoundedIntegerOption("number_look_ahead", "Sets limit of look-ahead strong-branching trials",
        0, 0,"");
    roptions->setOptionExtraInfo("number_look_ahead", 31);
  }


  void
  BonChooseVariable::computeMultipliers(double& upMult, double& downMult) const
  {
    const double* upTotalChange = pseudoCosts_.upTotalChange();
    const double* downTotalChange = pseudoCosts_.downTotalChange();
    const int* upNumber = pseudoCosts_.upNumber();
    const int* downNumber = pseudoCosts_.downNumber();
    double sumUp=0.0;
    double numberUp=0.0;
    double sumDown=0.0;
    double numberDown=0.0;
    for (int i=pseudoCosts_.numberObjects() - 1; i >= 0; --i) {
      sumUp += upTotalChange[i];
      numberUp += upNumber[i];
      sumDown += downTotalChange[i];
      numberDown += downNumber[i];
      message(PS_COST_HISTORY)
      <<i<< upNumber[i]<< upTotalChange[i]
      << downNumber[i]<< downTotalChange[i]<<CoinMessageEol;
    }
    upMult=(1.0+sumUp)/(1.0+numberUp);
    downMult=(1.0+sumDown)/(1.0+numberDown);

    message(PS_COST_MULT)
    <<upMult<< downMult<<CoinMessageEol;
  }

  double
  BonChooseVariable::computeUsefulness(const double MAXMIN_CRITERION,
      const double upMult, const double downMult,
      const double value,
      const OsiObject* object, int i,
      double& value2) const
  {
#ifdef OLD_USEFULLNESS
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
#else
   // FIXME: Hanlding initialization correctly
    int numberUp = pseudoCosts_.upNumber()[i];
    int numberDown = pseudoCosts_.downNumber()[i];
    if (sortCrit_ >= DecrPs && sortCrit_ <= IncrPs) {//Using pseudo costs
      double sumUp = pseudoCosts_.upTotalChange()[i]+1.0e-30;
      double sumDown = pseudoCosts_.downTotalChange()[i]+1.0e-30;
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
      if (numberUp == 0 || numberDown == 0) {
        if (value == 0.) useful = 0;
        else if (MAXMIN_CRITERION >= 0.5)//Sort on max infeasibility
          useful =  CoinMin(upEst, downEst);
        else {//Do min infeasibility
          useful = CoinMax(upEst, downEst);
        }
      }
      value2 = -COIN_DBL_MAX;
      if (numberUp   < numberBeforeTrustedList_ ||
          numberDown < numberBeforeTrustedList_) {
        value2 = value;
      }
#endif
      message(PS_COST_ESTIMATES)
      <<i<< useful<< upEst<< downEst<< value<< value2<< CoinMessageEol;
      return useful;
    }
#ifndef OLD_USEFULLNESS
    else if (sortCrit_ >= DecrInfeas && sortCrit_ <= IncrInfeas) {//Just return infeasibility
      double usefull = value;
      value2 = value;
      return usefull;
    }
    else {
      throw CoinError("BonChooseVariable", "computeUsefullnee",
          "Unknown criterion for soring candidates");
      return COIN_DBL_MAX;
    }
  }
#endif

  int
  BonChooseVariable::setupList ( OsiBranchingInformation *info, bool initialize)
  {
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
      throw CoinError(CNAME,"setupList","Should not be called with initialize==false");
    }
    numberOnList_=0;
    numberUnsatisfied_=0;
    int numberObjects = solver_->numberObjects();
    assert (numberObjects);
    if (numberObjects>pseudoCosts_.numberObjects()) {
      //std::cout<<"Number objects "<<numberObjects<<std::endl;
      //AW : How could that ever happen?
      //PB : It happens for instance when SOS constraints are added. They are added after the creation of this.
      //   assert(false && "Right now, all old content is deleted!");
      // redo useful arrays
      int saveNumberBeforeTrusted = pseudoCosts_.numberBeforeTrusted();
      pseudoCosts_.initialize(numberObjects);
      pseudoCosts_.setNumberBeforeTrusted(saveNumberBeforeTrusted);
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
              }
              else {
                check2=0.0;
                checkIndex2 = j;
                break;
              }
            }
          }
        }
        else {
          // worse priority
          // to end
          assert (list_[putOther-1]<0);
          list_[--putOther]=i;
          maximumStrong = CoinMin(maximumStrong,putOther);
        }
      }
    }
#if 0
    for (int i=0; i<maximumStrong; i++) {
      int way;
      message(CANDIDATE_LIST)<<i
      <<list_[i]<<i<<useful_[i]<<object[list_[i]]->infeasibility(info,way)
      <<CoinMessageEol;
    }
#endif
    // Get list
    numberOnList_=0;
    if (feasible) {
      maximumStrong = CoinMin(maximumStrong,putOther);
      for (i=0;i<maximumStrong;i++) {
      if (list_[i]>=0) {
#ifdef OLD_USEFULLNESS
	list_[numberOnList_]=list_[i];
	useful_[numberOnList_++]=-useful_[i];
      
#else
          list_[numberOnList_]=list_[i];
          if ((sortCrit_ & 1) == 0) {
            useful_[numberOnList_++]=-useful_[i];
          }
          else useful_[numberOnList_++] = useful_[i];
#endif
          message(CANDIDATE_LIST2)<<numberOnList_-1
          <<list_[numberOnList_-1]<<numberOnList_-1<<useful_[numberOnList_-1]
          <<CoinMessageEol;
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
              message(CANDIDATE_LIST3)<<number_not_trusted_-1
              <<list2[number_not_trusted_-1]<<number_not_trusted_-1
              <<useful2[number_not_trusted_-1]<<CoinMessageEol;
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
    }
    else {
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
      for (int i=0; i<numberOnList_; i++){
        message(CANDIDATE_LIST)<<i<< list_[i]<< i<< useful_[i]
        <<object[list_[i]]->infeasibility(info,way)
        <<CoinMessageEol;
        }
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
    bool isRoot = isRootNode(info);
    int numberStrong = numberStrong_;
    if (isRoot) {
      numberStrong = CoinMax(numberStrong_, numberStrongRoot_);
    }
    std::vector<double> save_sol;
    if (bb_log_level_>=3) {
       save_sol.resize(info->numberColumns_);
       std::copy(info->solution_, info->solution_ + info->numberColumns_ , save_sol.begin());
    }
    if (numberUnsatisfied_) {
      const double* upTotalChange = pseudoCosts_.upTotalChange();
      const double* downTotalChange = pseudoCosts_.downTotalChange();
      const int* upNumber = pseudoCosts_.upNumber();
      const int* downNumber = pseudoCosts_.downNumber();
      int numberBeforeTrusted = pseudoCosts_.numberBeforeTrusted();
      int numberLeft = CoinMin(numberStrong -numberStrongDone_,numberUnsatisfied_);
      results_.clear();
      int returnCode=0;
      bestObjectIndex_ = -1;
      bestWhichWay_ = -1;
      firstForcedObjectIndex_ = -1;
      firstForcedWhichWay_ =-1;
      double bestTrusted=-COIN_DBL_MAX;
      for (int i=0;i<numberLeft;i++) {
        int iObject = list_[i];
        if (numberBeforeTrusted == 0||
            i < minNumberStrongBranch_ ||
            (
              only_pseudo_when_trusted_ && number_not_trusted_>0 ) ||
              ( !isRoot && (upNumber[iObject]<numberBeforeTrusted ||
                          downNumber[iObject]<numberBeforeTrusted ))||
              ( isRoot && (!upNumber[iObject] && !downNumber[iObject])) ) {
         
             results_.push_back(HotInfo(solver, info,
                                solver->objects(), iObject));
        }
        else {
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
      if (results_.size() > 0) {
        returnCode = doStrongBranching(solver, info, (int)results_.size(), 1);
        if (bb_log_level_>=3) {
          OsiObject ** obj = solver->objects();
          const char* stat_msg[] = {"NOTDON", "FEAS", "INFEAS", "NOFINI"};
          message(SB_START)<<info->objectiveValue_<<CoinMessageEol;
          message(SB_HEADER)<<CoinMessageEol;
          for (unsigned int i = 0; i< results_.size(); i++) {
            double up_change = results_[i].upChange();
            double down_change = results_[i].downChange();
            int up_status = results_[i].upStatus();
            int down_status = results_[i].downStatus();
            int icol = obj[results_[i].whichObject()]->columnNumber();
            double val = save_sol[icol];
            message(SB_RES)<<(int) i<<icol<<val<<stat_msg[down_status+1]<<down_change
            <<stat_msg[up_status+1]<< up_change<< CoinMessageEol;
          }
        }
        if (returnCode>=0&&returnCode<=2) {
          if (returnCode) {
            returnCode=4;
            if (bestObjectIndex_>=0)
              returnCode=3;
          }
          for (unsigned int i=0;i < results_.size();i++) {
            int iObject = results_[i].whichObject();
            double upEstimate;
            if (results_[i].downStatus()== 2 || results_[i].upStatus()==2) {
              //continue;
            }
            if (results_[i].upStatus()!=1) {
              assert (results_[i].upStatus()>=0);
              upEstimate = results_[i].upChange();
            }
            else {
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
            }
            else {
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
        }
        else if (returnCode==3) {
          // max time - just choose one
          bestObjectIndex_ = list_[0];
          bestWhichWay_ = 0;
          returnCode=0;
        }
      }
      else {
        bestObjectIndex_=list_[0];
      }
      if ( bestObjectIndex_ >=0 ) {
        OsiObject * obj = solver->objects()[bestObjectIndex_];
        obj->setWhichWay(bestWhichWay_);
        message(BRANCH_VAR)<<obj->columnNumber()<< bestWhichWay_
        <<CoinMessageEol;
      }
      message(CHOSEN_VAR)<<bestObjectIndex_<<CoinMessageEol;
      if (numberFixed==numberUnsatisfied_&&numberFixed)
        returnCode=4;
      return returnCode;
    }
    else {
      return 1;
    }
  }

  /*  This is a utility function which does strong branching on
      a list of objects and stores the results in OsiHotInfo.objects.
      On entry the object sequence is stored in the OsiHotInfo object
      and maybe more.
      It returns -
      -1 - one branch was infeasible both ways
       0 - all inspected - nothing can be fixed
       1 - all inspected - some can be fixed (returnCriterion==0)
       2 - may be returning early - one can be fixed (last one done) (returnCriterion==1) 
       3 - returning because max time
  */
  int 
  BonChooseVariable::doStrongBranching( OsiSolverInterface * solver, 
  				    OsiBranchingInformation *info,
  				    int numberToDo, int returnCriterion)
  {
    // Prepare stuff for look-ahead heuristic
    double bestLookAhead_ = -COIN_DBL_MAX;
    int trialsSinceBest_ = 0;
    bool isRoot = isRootNode(info);
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
      if (status0==3) {
        // new solution already saved
        if (trustStrongForSolution_) {
        info->cutoff_ = goodObjectiveValue_;
        status0=0;
        }
      }
      if(solver->getRowCutDebugger() && status0 == 1 ){
           OsiTMINLPInterface * tminlp_solver = dynamic_cast<OsiTMINLPInterface *> (solver);
           throw tminlp_solver->newUnsolvedError(1, tminlp_solver->problem(), "SB");
      }
      numberStrongIterations_ += thisSolver->getIterationCount();
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
        fflush(stdout);
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
      if (status1==3) {
        // new solution already saved
        if (trustStrongForSolution_) {
        info->cutoff_ = goodObjectiveValue_;
        status1=0;
        }
      }
      if(solver->getRowCutDebugger() && status1 == 1){
           OsiTMINLPInterface * tminlp_solver = dynamic_cast<OsiTMINLPInterface *> (solver);
           throw tminlp_solver->newUnsolvedError(1, tminlp_solver->problem(), "SB");
      }
      numberStrongIterations_ += thisSolver->getIterationCount();
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
      if (status0==1&&status1==1) {
        // infeasible
        returnCode=-1;
        //break; // exit loop
      } else if (status0==1||status1==1) {
        numberStrongFixed_++;
  	returnCode=1;
      }
      bool hitMaxTime = ( CoinCpuTime()-timeStart > info->timeRemaining_)
                        || ( CoinCpuTime() - start_time_ > time_limit_);
      if (hitMaxTime) {
        returnCode=3;
        break;
      }
      // stop if look ahead heuristic tells us so
      if (!isRoot && numberLookAhead_) {
	assert(status0==0 && status1==0);
	double upEstimate = result->upChange();
	double downEstimate = result->downChange();
	double MAXMIN_CRITERION = maxminCrit(info);
	double value = MAXMIN_CRITERION*CoinMin(upEstimate,downEstimate) + (1.0-MAXMIN_CRITERION)*CoinMax(upEstimate,downEstimate);
	if (value > bestLookAhead_) {
	  bestLookAhead_ = value;
	  trialsSinceBest_ = 0;
	}
	else {
	  trialsSinceBest_++;
	  if (trialsSinceBest_ >= numberLookAhead_) {
	    break;
	  }
	}
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

  bool BonChooseVariable::isRootNode(const OsiBranchingInformation *info) const
  {
    return info->depth_ == 0;
  }

  double
  BonChooseVariable::maxminCrit(const OsiBranchingInformation *info) const
  {
    double retval = maxmin_crit_no_sol_;
    if (cbc_model_) {
      // FIXME: should be replaced by info->stateOfSearch_
      const int stateOfSearch = cbc_model_->stateOfSearch();
      const int depth = info->depth_;
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
    if(!trustStrongForPseudoCosts_) return;
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
      //if (hotInfo->upStatus()!=1) 
      // AW: Let's update the pseudo costs only if the strong branching
      // problem was marked as "solved"
      if (hotInfo->upStatus()==0) {
        assert (hotInfo->upStatus()>=0);
        upTotalChange[index] += hotInfo->upChange()/object->upEstimate();
        upNumber[index]++;
      }
      else if (hotInfo->upStatus()==1) {
        // infeasible - just say expensive
        upNumber[index]++;
        if (info->cutoff_<1.0e50)
          upTotalChange[index] += 2.0*(info->cutoff_-info->objectiveValue_)/object->upEstimate();
        else
          upTotalChange[index] += 2.0*fabs(info->objectiveValue_)/object->upEstimate();
      }
    }
    else {
      if (hotInfo->downStatus()==0) {
        assert (hotInfo->downStatus()>=0);
        downTotalChange[index] += hotInfo->downChange()/object->downEstimate();
        downNumber[index]++;
      }
      else if (hotInfo->downStatus()==1) {
        downNumber[index]++;
        // infeasible - just say expensive
        if (info->cutoff_<1.0e50)
          downTotalChange[index] += 2.0*(info->cutoff_-info->objectiveValue_)/object->downEstimate();
        else
          downTotalChange[index] += 2.0*fabs(info->objectiveValue_)/object->downEstimate();
      }
    }
  }

// Given a branch fill in useful information e.g. estimates 
void  
BonChooseVariable::updateInformation( int index, int branch,  
                                      double changeInObjective, double changeInValue, 
                                      int status) 
{ 
  if(cbc_model_ == NULL) return;
  assert (index<solver_->numberObjects()); 
  assert (branch<2); 

  if(fabs(changeInValue) < 1e-6) return;

  double* upTotalChange = pseudoCosts_.upTotalChange(); 
  double* downTotalChange = pseudoCosts_.downTotalChange(); 
  int* upNumber = pseudoCosts_.upNumber(); 
  int* downNumber = pseudoCosts_.downNumber(); 
    message(UPDATE_PS_COST)<<index<< branch
    <<changeInObjective<<changeInValue<<status
    <<CoinMessageEol;

  if (branch) { 
    if (status!=1) { 
      assert (status>=0); 
      upTotalChange[index] += changeInObjective/changeInValue; 
      upNumber[index]++; 
    } else { 
      // infeasible - just say expensive 
      upNumber[index]++; 
      assert(cbc_model_); // Later, we need to get this information in a different way... 
      double cutoff = cbc_model_->getCutoff(); 
      double objectiveValue = cbc_model_->getCurrentObjValue(); 
      if (cutoff<1.0e50) 
        upTotalChange[index] += 2.0*(cutoff-objectiveValue)/changeInValue; 
      else 
        upTotalChange[index] += 2.0*fabs(objectiveValue)/changeInValue; 
    } 
  } else { 
    if (status!=1) { 
      assert (status>=0); 
      downTotalChange[index] += changeInObjective/changeInValue; 
      downNumber[index]++; 
    } else { 
      assert(cbc_model_); 
      // infeasible - just say expensive 
      downNumber[index]++; 
      double cutoff = cbc_model_->getCutoff(); 
      double objectiveValue = cbc_model_->getCurrentObjValue(); 
      if (cutoff<1.0e50) 
        downTotalChange[index] += 2.0*(cutoff-objectiveValue)/changeInValue; 
      else 
        downTotalChange[index] += 2.0*fabs(objectiveValue)/changeInValue; 
    } 
  }   
} 


  HotInfo::HotInfo(): OsiHotInfo(), infeasibilities_(){
  } 

  HotInfo::HotInfo( OsiSolverInterface * solver,
                    const OsiBranchingInformation *info,
                    const OsiObject * const * objects,
                    int whichObject): 
  OsiHotInfo(solver, info, objects, whichObject),
  infeasibilities_(){
     infeasibilities_.resize(branchingObject_->numberBranches());
  }

  HotInfo::HotInfo(const HotInfo& other): OsiHotInfo(other),
                                          infeasibilities_(other.infeasibilities_){
  }

  OsiHotInfo * 
  HotInfo::clone() const {
    return new HotInfo(*this);
  }

  HotInfo &
  HotInfo::operator=(const HotInfo &rhs){
    if(this != &rhs){
      OsiHotInfo::operator=(rhs);
      infeasibilities_ = rhs.infeasibilities_;
    }
    return (*this);
  }

  HotInfo::~HotInfo(){
  }

  int
  HotInfo::updateInformation(const OsiSolverInterface * solver,
                             const OsiBranchingInformation * info,
                             OsiChooseVariable * choose){
    //printf("in HotInfo::updateInformation\n");
    int iBranch = branchingObject_->branchIndex()-1;
    double & infeasibility = infeasibilities_[iBranch] = 0.;
   
    OsiObject ** objects = solver->objects();
    int numObject = solver->numberObjects();
    for(int i = 0 ; i < numObject ; i++){
       infeasibility += objects[i]->checkInfeasibility(info);
    }
    int status = OsiHotInfo::updateInformation(solver, info, choose);
#if 1
    if(!solver->isProvenPrimalInfeasible() && !solver->isProvenOptimal()){
      status = 2;
      statuses_[iBranch] = status;
    }
    else if(solver->isProvenPrimalInfeasible() && fabs(solver->getObjValue()) < 1e-06) {
      assert(solver->messageHandler() != NULL);
      *solver->messageHandler() << "Very small infeasibility: " << solver->getObjValue() << CoinMessageEol;
      status = 2;
      statuses_[iBranch] = status;
    }
#endif
    return status;
  }
  
}/* Ends Bonmin's namespace.*/


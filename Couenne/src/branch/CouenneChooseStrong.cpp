/*
 * Name:    CouenneChooseStrong.cpp
 * Authors: Andreas Waechter, IBM Corp.
 * Purpose: Strong branching objects for Couenne
 *
 * (C) Carnegie-Mellon University, 2006-08.
 * This file is licensed under the Common Public License (CPL)
 */

#include "CoinTime.hpp"
#include "CouenneChooseStrong.hpp"
#include "CouenneProblem.hpp"
#include "CouenneBranchingObject.hpp"

namespace Bonmin {

  /// constructor
  CouenneChooseStrong::CouenneChooseStrong (BabSetupBase &b, CouenneProblem* p, JnlstPtr jnlst) :

    BonChooseVariable (b, b.continuousSolver()),
    problem_          (p),
    jnlst_            (jnlst) {

    std::string s;
    b.options () -> GetStringValue ("pseudocost_mult_lp", s, "couenne.");
    pseudoUpdateLP_ = (s == "yes");      
  }

  /// copy constructor
  CouenneChooseStrong::CouenneChooseStrong (const CouenneChooseStrong& rhs) :
    BonChooseVariable (rhs),
    problem_          (rhs.problem_),
    pseudoUpdateLP_   (rhs.pseudoUpdateLP_),
    jnlst_            (rhs.jnlst_)
  {}

  /// destructor
  CouenneChooseStrong::~CouenneChooseStrong()
  {}

  /// cloning method
  OsiChooseVariable *
  CouenneChooseStrong::clone() const
  {
    return new CouenneChooseStrong(*this);
  }

  /// assignment operator
  CouenneChooseStrong&
  CouenneChooseStrong::operator=(const CouenneChooseStrong & rhs)
  {
    if (this != &rhs) {
      BonChooseVariable::operator=(rhs);
      problem_ = rhs.problem_;
    }
    return *this;
  }

  /// compute Euclidean distance between two points (most likely LP solutions)
  /// l_2 norm by default, but can change it by fourth parameter
  double distance (const double *p1, const double *p2, register int size, double k=2.) {

    register double 
      result = 0,
      element;

    while (size--) {
      element = *p1++ - *p2++;
      result += pow (element, k);
    }

    return pow (result, 1./k);
  }

  /**  This is a utility function which does strong branching on
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
  int CouenneChooseStrong::doStrongBranching (OsiSolverInterface * solver, 
					      OsiBranchingInformation *info,
					      int numberToDo, int returnCriterion)
  {

    jnlst_ -> Printf (J_DETAILED, J_BRANCHING, "\n-\n------- CCS: trying %d objects:\n", numberToDo);

    int numberColumns = solver -> getNumCols ();

    solver -> markHotStart (); // save current LP point

    const double
      *lower = info -> lower_,
      *upper = info -> upper_;

    double 
      *saveLower = CoinCopyOfArray (info -> lower_, numberColumns),
      *saveUpper = CoinCopyOfArray (info -> upper_, numberColumns),
      *Lower0    = CoinCopyOfArray (info -> lower_, numberColumns), // delete afterwards
      *Upper0    = CoinCopyOfArray (info -> upper_, numberColumns),
      *oldLower  = new double [numberColumns],
      *oldUpper  = new double [numberColumns],
      *lpSol     = NULL, 
       timeStart = CoinCpuTime ();

    // LP solution for distance
    if (pseudoUpdateLP_) 
      lpSol = CoinCopyOfArray (info -> solution_, numberColumns);

    // provide Couenne problem with point/bounds contained in info
    problem_ -> domain () -> push
      (problem_ -> nVars (),
       info -> solution_,
       info -> lower_,
       info -> upper_);

    int returnCode = 0, iDo = 0;

    for (;iDo < numberToDo; iDo++) {

      HotInfo * result = results_ () + iDo; // retrieve i-th object to test

      CouenneObject *CouObj = dynamic_cast <CouenneObject *>
	(solver_ -> objects () [result -> whichObject ()]);

      // For now just 2 way
      OsiBranchingObject * branch = result -> branchingObject ();
      assert (branch->numberBranches()==2);

      CouenneBranchingObject *cb = dynamic_cast <CouenneBranchingObject *> (branch);
      if (cb) cb -> setSimulate (true);

      /* Try the first direction.  Each subsequent call to branch()
	 performs the specified branch and advances the branch object
	 state to the next branch alternative.) */

      int 
	status0 = -1, 
	status1 = -1;

      OsiSolverInterface * thisSolver = solver; 

      // DOWN DIRECTION ///////////////////////////////////////////////////////

      if (branch -> boundBranch ()) { // a (variable) bound branch

        if (branch -> branch (solver) > COUENNE_INFINITY) // branch is infeasible
	  result -> setDownStatus (status0 = 1);

	else { // branch is feasible, solve and compare

	  solver -> solveFromHotStart ();

	  if (pseudoUpdateLP_ && CouObj)
	    CouObj -> setEstimate (distance (lpSol, solver -> getColSolution (), numberColumns), 0);
	}

      } else {                       // some more complex branch, have to clone solver

        // adding cuts or something 
        thisSolver = solver -> clone ();
        if (branch -> branch (thisSolver) > COUENNE_INFINITY)
	  result -> setDownStatus (status0 = 1);

	else { // set hot start iterations
	  int limit;
	  thisSolver -> getIntParam (OsiMaxNumIterationHotStart, limit);
	  thisSolver -> setIntParam (OsiMaxNumIteration,         limit); 

	  thisSolver -> resolve ();
	  if (pseudoUpdateLP_ && CouObj)
	    CouObj -> setEstimate (distance (lpSol, thisSolver->getColSolution (), numberColumns), 0);
	}
      }

      // can check if we got solution
      // status is 0 finished, 1 infeasible and 2 unfinished and 3 is solution

      if (status0 < 0) {
	status0 = result->updateInformation(thisSolver,info,this);
	numberStrongIterations_ += thisSolver->getIterationCount();
      }

      if ((status0 == 3) && (trustStrongForSolution_)) {
        // new solution already saved
  	info -> cutoff_ = goodObjectiveValue_;
	problem_ -> setCutOff (goodObjectiveValue_);
  	status0 = 0;
      }

      if (solver != thisSolver)
        delete thisSolver;

      // save current bounds as tightened by the down branch; will be
      // used below to update global bounding box in solver
      CoinCopyN (problem_ -> Lb (), numberColumns, oldLower);
      CoinCopyN (problem_ -> Ub (), numberColumns, oldUpper);

      // Restore pre-left-branch bounds in solver
      for (int j=0; j<numberColumns; j++) {

        if (saveLower [j] != lower [j]) solver -> setColLower (j, saveLower [j]);
        if (saveUpper [j] != upper [j]) solver -> setColUpper (j, saveUpper [j]);
      }

      // UP DIRECTION ///////////////////////////////////////////////////////

      thisSolver = solver; 

      if (branch -> boundBranch ()) { // (variable) bound branch 

        if (branch -> branch (solver) > COUENNE_INFINITY)
	  result -> setUpStatus (status1 = 1);

        else {
	  solver -> solveFromHotStart ();
	  if (pseudoUpdateLP_ && CouObj) 
	    CouObj -> setEstimate (distance (lpSol, solver -> getColSolution (), numberColumns), 1);
	}
      } else {                     // some more complex branch, have to clone solver
        // adding cuts or something 
        thisSolver = solver -> clone ();
        if (branch -> branch (thisSolver) > COUENNE_INFINITY)
	  result -> setUpStatus (status1 = 1);

	else {
	  // set hot start iterations
	  int limit;
	  thisSolver -> getIntParam (OsiMaxNumIterationHotStart, limit);
	  thisSolver -> setIntParam (OsiMaxNumIteration,         limit); 

	  thisSolver -> resolve();
	  if (pseudoUpdateLP_ && CouObj) 
	    CouObj -> setEstimate (distance (lpSol, thisSolver->getColSolution (), numberColumns), 1);
	}
      }

      // can check if we got solution
      // status is 0 finished, 1 infeasible and 2 unfinished and 3 is solution

      if (status1 < 0) {
	status1 = result->updateInformation(thisSolver,info,this);
	numberStrongDone_++;
      }

      numberStrongIterations_ += thisSolver->getIterationCount();

      if ((status1==3) && (trustStrongForSolution_)) {
        // new solution already saved
	info -> cutoff_ = goodObjectiveValue_;
	problem_ -> setCutOff (goodObjectiveValue_);
	status1 = 0;
      }

      if (cb) cb -> setSimulate (false);

      /////////////////////////////////////////////////////////////////////////////

      if (solver != thisSolver)
        delete thisSolver;

      bool tightened = false;

      t_chg_bounds *chg_bds = new t_chg_bounds [numberColumns];

      // extend problem_'s bounding box to include downbranch's tightened
      for (int j=0; j<numberColumns; j++) {

        if (oldLower [j] < problem_ -> Lb (j)) problem_ -> Lb (j) = oldLower [j];
        if (oldUpper [j] > problem_ -> Ub (j)) problem_ -> Ub (j) = oldUpper [j];

	if (problem_ -> Lb (j) > lower [j] + COUENNE_EPS) {
	  chg_bds [j].setLower (t_chg_bounds::CHANGED);
	  tightened = true;
	}

	if (problem_ -> Ub (j) < upper [j] - COUENNE_EPS) {
	  chg_bds [j].setUpper (t_chg_bounds::CHANGED);
	  tightened = true;
	}
      }

      if (tightened &&                     // have tighter bounds
	  !(problem_ -> btCore (chg_bds))) // tighten again on root

	status0 = status1 = 1;	           // if returns false, problem is infeasible

      delete [] chg_bds;

      // create union of bounding box from both branching directions
      for (int j=0; j<numberColumns; j++) {

        if (oldLower [j] < problem_ -> Lb (j)) problem_ -> Lb (j) = oldLower [j];
        if (oldUpper [j] > problem_ -> Ub (j)) problem_ -> Ub (j) = oldUpper [j];
      }

      // set new bounding box as the possibly tightened one (a subset
      // of the initial)
      for (int j=0; j<numberColumns; j++) {

        solver -> setColLower (j, saveLower [j] = problem_ -> Lb (j));
        solver -> setColUpper (j, saveUpper [j] = problem_ -> Ub (j));
      }

      /*
        End of evaluation for this candidate object. Possibilities are:

        * Both sides below cutoff; this variable is a candidate for
          branching.

        * Both sides infeasible or above the objective cutoff: no
          further action here. Break from the evaluation loop and
          assume the node will be purged by the caller.

        * One side feasible and below cutoff: Install the branch
          (i.e., fix the variable). Possibly break from the evaluation
          loop and assume the node will be reoptimised by the caller.
      */

      if (status0 == 1 && 
	  status1 == 1) { // infeasible
        returnCode=-1;
        break; // exit loop
      } else if (status0==1 || status1==1) {
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
    } // end loop /***********************************/
  

    if (jnlst_ -> ProduceOutput (J_DETAILED, J_BRANCHING)) {
      printf ("tightened bounds: ");
      // create union of bounding box from both branching directions
      for (int j=0; j<numberColumns; j++) {
      
	if (problem_ -> Lb (j) > Lower0 [j]) printf ("l%d (%g-->%g) ", j,Lower0[j], problem_->Lb (j));
	if (problem_ -> Ub (j) < Upper0 [j]) printf ("u%d (%g-->%g) ", j,Upper0[j], problem_->Ub (j));
      }
    }

    delete [] Lower0;
    delete [] Upper0;

    problem_ -> domain () -> pop (); // discard current point/bounds from problem

    delete [] lpSol;

    jnlst_ -> Printf (J_DETAILED, J_BRANCHING, "----------------------done\n\n\n");

    if (iDo < numberToDo) iDo++; // exited due to infeasibility
    assert (iDo <= (int) results_.size());
    results_.resize (iDo);

    delete [] saveLower;
    delete [] saveUpper;

    solver -> unmarkHotStart ();     // Delete the snapshot

    return returnCode;
  }


#ifndef CCS_EXPERIMENTAL
  // Sets up strong list and clears all if initialize is true.
  // Returns number of infeasibilities.
  int CouenneChooseStrong::setupList (OsiBranchingInformation *info, bool initialize) {

    initialize = true; // to avoid failed assert in BonChooseVariable::setupList()

    problem_ -> domain () -> push 
      (problem_ -> nVars (),
       info -> solution_, 
       info -> lower_, 
       info -> upper_); // have to alloc+copy

    if (jnlst_ -> ProduceOutput (J_DETAILED, J_BRANCHING)) {
      printf ("----------------- (strong) setup list\n");
      for (int i=0; i<problem_ -> domain () -> current () -> Dimension (); i++)
	printf ("%4d %20.4g [%20.4g %20.4g]\n", i,
		info -> solution_ [i], info -> lower_ [i], info -> upper_ [i]);
    }

    // call Bonmin's setuplist
    int retval = BonChooseVariable::setupList (info, initialize);

    jnlst_ -> Printf (J_DETAILED, J_BRANCHING, "----------------- (strong) setup list done\n");

    problem_ -> domain () -> pop ();
    return retval;
  }
#endif


  /// Add list of options to be read from file ////////////////////////////////////////
  void CouenneChooseStrong::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

    roptions -> AddStringOption6
      ("pseudocost_mult",
       "Multipliers of pseudocosts for estimating and update estimation of bound",
       "infeasibility",

       "infeasibility", "infeasibility returned by object",

       "projectDist",   "distance between current LP point and resulting branches' LP points",

       "interval_lp",   "width of the interval between bound and current lp point",
       "interval_lp_rev",   "similar to interval_lp, reversed",

       "interval_br",   "width of the interval between bound and branching point",
       "interval_br_rev",   "similar to interval_br, reversed");

    roptions -> AddStringOption2
      ("pseudocost_mult_lp",
       "Use distance between LP points to update multipliers of pseudocosts "  
       "after simulating branching",
       "no",
       "yes", "",
       "no",  "");
  }


  // Returns true if solution looks feasible against given objects
  bool CouenneChooseStrong::feasibleSolution (const OsiBranchingInformation * info,
					      const double * solution,
					      int numberObjects,
					      const OsiObject ** objects)

  {return problem_ -> checkNLP (solution, solution [problem_ -> Obj (0) -> Body () -> Index ()]);}


#ifdef CCS_EXPERIMENTAL

  ///
  int CouenneChooseStrong::setupList (OsiBranchingInformation *info, bool initialize) {

    if (jnlst_ -> ProduceOutput (J_DETAILED, J_BRANCHING)) {
      printf ("----------------- (strong) setup list on %d objects\n",
	      solver_ -> numberObjects());
      for (int i=0; i<problem_ -> domain () -> current () -> Dimension (); i++)
	printf ("%4d %20.4g [%20.4g %20.4g]\n", i,
		info -> solution_ [i], info -> lower_ [i], info -> upper_ [i]);
    }

    if (numberBeforeTrustedList_ < 0) {
      number_not_trusted_ = 1;
      return OsiChooseVariable::setupList(info, initialize);
    }

    initialize = true; // to avoid failed assert in BonChooseVariable::setupList()

    problem_ -> domain () -> push 
      (problem_ -> nVars (),
       info -> solution_, 
       info -> lower_, 
       info -> upper_); // have to alloc+copy

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

    numberOnList_=0;
    numberUnsatisfied_=0;
    int numberObjects = solver_->numberObjects();
    assert (numberObjects);
    if (numberObjects>pseudoCosts_.numberObjects()) {
      //std::cout<<"Number objects "<<numberObjects<<std::endl;
      //AW : How could that ever happen?  Right now, all old content is deleted!
      //   assert(false && "Right now, all old content is deleted!");
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
    int max_most_fra = setup_pseudo_frac_ > 0. ? 
      (int)floor(setup_pseudo_frac_*(double)maximumStrong): 0;

    if (setup_pseudo_frac_ > 0.)
      max_most_fra = CoinMax(1, max_most_fra);

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

    // 1. check the infeasibility of all infeas objects
    // 2. use match infeas <-> branch to assign infeasibility to branch objects

    /// scan all branch objects
    for ( i=0;i<numberObjects;i++) {
      int way;

      double value = object [i] -> infeasibility (info, way); // branch object's infeas

      if (value>0.0) {

        numberUnsatisfied_++;

        if (value>=1e50) { // should never be
          // infeasible
          feasible=false;
          break;
        }

        int priorityLevel = object [i]->priority();

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

    // Get list
    numberOnList_=0;
    if (feasible) {
      maximumStrong = CoinMin(maximumStrong,putOther);
      for (i=0;i<maximumStrong;i++) {
        if (list_[i]>=0) {
          list_[numberOnList_]=list_[i];
          if ((sortCrit_ & 1) == 0) {
            useful_[numberOnList_++]=-useful_[i];
          }
          else useful_[numberOnList_++] = useful_[i];
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
      for (int i=0; i<numberOnList_; i++)
        message(CANDIDATE_LIST)<<i<< list_[i]<< i<< useful_[i]
        <<object[list_[i]]->infeasibility(info,way)
        <<CoinMessageEol;
    }

    jnlst_ -> Printf (J_DETAILED, J_BRANCHING, "----------------- (strong) setup list done\n");

    problem_ -> domain () -> pop ();

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
  int CouenneChooseStrong::chooseVariable (OsiSolverInterface * solver,
					   OsiBranchingInformation *info,
					   bool fixVariables) {

    // We assume here that chooseVariable is called once at the very
    // beginning with fixVariables set to true.  This is then the root
    // node.
    bool isRoot = isRootNode(info);
    int numberStrong = numberStrong_;
    if (isRoot) {
      numberStrong = CoinMax(numberStrong_, numberStrongRoot_);
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
              !isRoot && (upNumber[iObject]<numberBeforeTrusted ||
                          downNumber[iObject]<numberBeforeTrusted )||
              isRoot && (!upNumber[iObject] && !downNumber[iObject]) ) {
         
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
        returnCode = doStrongBranching(solver, info, results_.size(), 1);
        if (bb_log_level_>=3) {
          const char* stat_msg[] = {"NOTDON", "FEAS", "INFEAS", "NOFINI"};
          message(SB_HEADER)<<CoinMessageEol;
          for (unsigned int i = 0; i< results_.size(); i++) {
            double up_change = results_[i].upChange();
            double down_change = results_[i].downChange();
            int up_status = results_[i].upStatus();
            int down_status = results_[i].downStatus();
            message(SB_RES)<<(int) i<<stat_msg[down_status+1]<<down_change
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
        obj->setWhichWay(	bestWhichWay_);
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


  // Given a candidate  fill in useful information e.g. estimates
  void CouenneChooseStrong::updateInformation (const OsiBranchingInformation *info,
					       int branch, OsiHotInfo * hotInfo) {

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
      else if (hotInfo->upStatus()==1) {
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
  void CouenneChooseStrong::updateInformation (int index, int branch,  
					       double changeInObjective, double changeInValue, 
					       int status) { 
    if(cbc_model_ == NULL) return;
    assert (index<solver_->numberObjects()); 
    assert (branch<2); 
    assert (changeInValue>0.0); 
    assert (branch<2); 
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
#endif
}

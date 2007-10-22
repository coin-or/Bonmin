// (C) Copyright International Business Machines Corporation 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Laszlo Ladanyi, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University

#include "CoinHelperFunctions.hpp"
#include "BCP_lp_node.hpp"
#include "BCP_lp.hpp"
#include "BM.hpp"
#include "BonChooseVariable.hpp"
#include "BonCurvBranchingSolver.hpp"
#include "BonQpBranchingSolver.hpp"
#include "BonLpBranchingSolver.hpp"
#include "BonOsiTMINLPInterface.hpp"

static bool ifprint = true;
static bool ifprint2 = false;
static bool ifprint3 = false;
// #define BM_PRINT_DATA

//#############################################################################

BCP_branching_decision
BM_lp::select_branching_candidates(const BCP_lp_result& lpres,
                                   const BCP_vec<BCP_var*>& vars,
                                   const BCP_vec<BCP_cut*>& cuts,
                                   const BCP_lp_var_pool& local_var_pool,
                                   const BCP_lp_cut_pool& local_cut_pool,
                                   BCP_vec<BCP_lp_branching_object*>& cands,
				   bool force_branch)
{
  /* FIXME: this is good only for BB and using NLP to solve a node */
  bm_stats.incNumberNodeSolves();

  Bonmin::OsiTMINLPInterface* nlp =
    dynamic_cast<Bonmin::OsiTMINLPInterface*>(getLpProblemPointer()->lp_solver);
  // If we are doing pure B&B then we have an nlp, and then we check for
  // consecutive failures.
  if (nlp) {
    if (lpres.termcode() & BCP_Abandoned) {
      if (nlp->isIterationLimitReached()) {
	print(ifprint, "\
BM_lp: At node %i : WARNING: nlp reached iter limit. Will force branching\n",
	       current_index());
      } else {
	print(ifprint, "\
BM_lp: At node %i : WARNING: nlp is abandoned. Will force branching\n",
	       current_index());
      }
      // nlp failed
      nlp->forceBranchable();
      numNlpFailed_++;
      if (numNlpFailed_ >= par.entry(BM_par::NumNlpFailureMax)) {
	print(ifprint, "WARNING! Too many (%i) NLP failures in a row. Abandoning node.",
	       numNlpFailed_);
	return BCP_DoNotBranch_Fathomed;
      }
    } else {
      numNlpFailed_ = 0;
    }
  }

  OsiBranchingInformation brInfo(nlp, false, true);
  brInfo.cutoff_ = upper_bound() + get_param(BCP_lp_par::Granularity);
  brInfo.objectiveValue_ = lpres.objval();
  brInfo.integerTolerance_ = integerTolerance_;
  brInfo.timeRemaining_ = get_param(BCP_lp_par::MaxRunTime) - CoinCpuTime();
  brInfo.numberSolutions_ = 0; /*FIXME*/
  brInfo.numberBranchingSolutions_ = 0; /*FIXME numBranchingSolutions_;*/
  brInfo.depth_ = current_level();

  BCP_branching_decision brDecision;
  if (bonmin_.getAlgorithm() == 0) {
    /* if pure B&B */
    brDecision = bbBranch(brInfo, cands);
  } else {
    brDecision = hybridBranch();
  }

  return brDecision;
}

//-----------------------------------------------------------------------------

BCP_branching_decision
BM_lp::hybridBranch()
{
    // FIXME: most of the pureBB stuff should work here.
    throw BCP_fatal_error("BM_lp: FIXME: make hybrid work...");
}

/*****************************************************************************/

void
BM_lp::unpack_pseudo_costs(BCP_buffer& buf)
{
  Bonmin::BonChooseVariable* choose =
    dynamic_cast<Bonmin::BonChooseVariable*>(bonmin_.branchingMethod());
  OsiPseudoCosts& pseudoCosts = choose->pseudoCosts();
  int numObj = pseudoCosts.numberObjects();
  double* upTotalChange = pseudoCosts.upTotalChange();
  int* upNumber = pseudoCosts.upNumber();
  double* downTotalChange = pseudoCosts.downTotalChange();
  int* downNumber = pseudoCosts.downNumber();
  
  buf.unpack(upTotalChange, numObj, false);
  buf.unpack(upNumber, numObj, false);
  buf.unpack(downTotalChange, numObj, false);
  buf.unpack(downNumber, numObj, false);

#ifdef BM_PRINT_DATA
  if (current_level() < 4) {
    print(ifprint2, "Received pseudocosts:\n");
    for (int i = 0; i < numObj; ++i) {
	const OsiObject* object = bonmin_.nonlinearSolver()->objects()[i];
	print(ifprint2, "col: %i,  obj: %i, dTot: %lf, dNum: %i, uTot: %lf, uNum: %i\n",
	      object->columnNumber(), i,
	      downTotalChange[i], downNumber[i],
	      upTotalChange[i], upNumber[i]);
    }
  }
#endif  
}


//-----------------------------------------------------------------------------

void
BM_lp::send_pseudo_cost_update(OsiBranchingInformation& branchInfo)
{
  bm_buf.clear();
  int itmp;
  double objchange;
  itmp = BM_PseudoCostUpdate;
  bm_buf.pack(itmp);
  for (int i = 0; i < objNum_; ++i) {
    const BM_SB_result& sbres = sbResult_[i];
    if ((sbres.branchEval & 1) != 0 && sbres.status[0] != BCP_Abandoned) {
      bm_buf.pack(sbres.objInd);
      itmp = 0;
      bm_buf.pack(itmp);
      if (sbres.status[0] == BCP_ProvenOptimal) {
	objchange = sbres.objval[0] - branchInfo.objectiveValue_;
      } else { // Must be BCP_ProvenPrimalInf
	if (branchInfo.cutoff_ < 1e50) {
	  objchange = 2*(branchInfo.cutoff_-branchInfo.objectiveValue_);
	} else {
	  objchange = 2*fabs(branchInfo.objectiveValue_);
	}
      }
      bm_buf.pack(objchange/sbres.varChange[0]);
    }
    if ((sbres.branchEval & 2) != 0 && sbres.status[1] != BCP_Abandoned) {
      bm_buf.pack(sbres.objInd);
      itmp = 1;
      bm_buf.pack(itmp);
      if (sbres.status[1] == BCP_ProvenOptimal) {
	objchange = sbres.objval[1] - branchInfo.objectiveValue_;
      } else { // Must be BCP_ProvenPrimalInf
	if (branchInfo.cutoff_ < 1e50) {
	  objchange = 2*(branchInfo.cutoff_-branchInfo.objectiveValue_);
	} else {
	  objchange = 2*fabs(branchInfo.objectiveValue_);
	}
      }
      bm_buf.pack(objchange/sbres.varChange[1]);
    }
  }
  itmp = -1;
  bm_buf.pack(itmp);
  send_message(parent(), bm_buf);

}

//-----------------------------------------------------------------------------

int
BM_lp::sort_objects(OsiBranchingInformation& branchInfo,
		    Bonmin::BonChooseVariable* choose, int& branchNum)
{
  const OsiObject* const * objects = branchInfo.solver_->objects();
  double upMult, downMult;
  choose->computeMultipliers(upMult, downMult);
  const double MAXMIN = choose->maxminCrit(&branchInfo);

  /* Order all objects that can be branched on */
  int lastPriority = objects[objInd_[0]]->priority();
  int infBlockStart = 0;
  int feasBlockStart = 0;
  branchNum = 0;
  infNum_ = 0;
  feasNum_ = 0;

  const bool isRoot = (current_index() == 0);
  int way;

  for (int i = 0; i < objNum_; ++i) {
    const int ind = objInd_[i];
    const OsiObject* object = objects[ind];
    double value = object->infeasibility(&branchInfo, way);
    if (value  > 0.0) {
      if (value >= 1e50) { // infeasible
	return -1;
      }
      int priorityLevel = object->priority();
      if (lastPriority < priorityLevel) {
	// sort the entries based on their usefulness
	if (infBlockStart < infNum_ &&
	    ! par.entry(BM_par::DisregardPriorities)) {
	  if (par.entry(BM_par::DecreasingSortInSetupList)) {
	    CoinSort_2(infUseful_ + infBlockStart, infUseful_ + infNum_,
		       infInd_ + infBlockStart,
		       CoinFirstGreater_2<double,int>());
	  } else {
	    CoinSort_2(infUseful_ + infBlockStart, infUseful_ + infNum_,
		       infInd_ + infBlockStart);
	  }
	}
	lastPriority = priorityLevel;
      }
      double dummy;
      infInd_[infNum_] = ind;
      if (par.entry(BM_par::UsePseudoCosts)) {
	infUseful_[infNum_] = isRoot ?
	  value : choose->computeUsefulness(MAXMIN, upMult, downMult, value,
					    object, ind, dummy);
      } else {
	infUseful_[infNum_] = value;
      }

      ++infNum_;
      branchNum += 2;

    } else { /* value == 0.0 */
      const OsiSOS* sos = dynamic_cast<const OsiSOS*>(object);
      if (sos) {
	// if an sos is feasible thne we don't do strong branching on it
	continue;
      }
      const int iCol = object->columnNumber();
      const double lb = branchInfo.lower_[iCol];
      const double ub = branchInfo.upper_[iCol];
      if (fabs(ub - lb) < 0.1) {
	continue;
      }
      value = branchInfo.solution_[iCol];
      ++branchNum;
      if (floor(value+0.5) > lb && ceil(value-0.5) < ub) {
	// The variable is integer, but neither at its lower nor at its upper
	// bound (the test accounts for tolerances)
	++branchNum;
      }
      int priorityLevel = object->priority();
      if (lastPriority < priorityLevel) {
	// sort the entries based on their usefulness
	if (feasBlockStart < feasNum_) {
	  if (par.entry(BM_par::DecreasingSortInSetupList) &&
	      ! par.entry(BM_par::DisregardPriorities)) {
	    CoinSort_2(feasUseful_ + feasBlockStart, feasUseful_ + feasNum_,
		       feasInd_ + feasBlockStart,
		       CoinFirstGreater_2<double,int>());
	  } else {
	    CoinSort_2(feasUseful_ + feasBlockStart, feasUseful_ + feasNum_,
		       feasInd_ + feasBlockStart);
	  }
	}
	lastPriority = priorityLevel;
      }
      double dummy;
      feasInd_[feasNum_] = ind;
      feasUseful_[feasNum_] = choose->computeUsefulness(MAXMIN,
							upMult, downMult, value,
							object, ind, dummy);
      ++feasNum_;
    }
  }
  if (infBlockStart < infNum_) {
    if (par.entry(BM_par::DecreasingSortInSetupList)) {
      CoinSort_2(infUseful_ + infBlockStart, infUseful_ + infNum_,
		 infInd_ + infBlockStart,
		 CoinFirstGreater_2<double,int>());
    } else {
      CoinSort_2(infUseful_ + infBlockStart, infUseful_ + infNum_,
		 infInd_ + infBlockStart);
    }
  }
  if (feasBlockStart < feasNum_) {
    if (par.entry(BM_par::DecreasingSortInSetupList)) {
      CoinSort_2(feasUseful_ + feasBlockStart, feasUseful_ + feasNum_,
		 feasInd_ + feasBlockStart,
		 CoinFirstGreater_2<double,int>());
    } else {
      CoinSort_2(feasUseful_ + feasBlockStart, feasUseful_ + feasNum_,
		 feasInd_ + feasBlockStart);
    }
  }
#ifdef BM_PRINT_DATA
  if (current_level() < 4) {
    OsiPseudoCosts& pseudoCosts = choose->pseudoCosts();
    double* upTotalChange = pseudoCosts.upTotalChange();
    int* upNumber = pseudoCosts.upNumber();
    double* downTotalChange = pseudoCosts.downTotalChange();
    int* downNumber = pseudoCosts.downNumber();
    
    print(ifprint2, "Before SB Infeas\n");
    for (int i = 0; i < infNum_; ++i) {
      const int ind = infInd_[i];
      const OsiObject* object = objects[ind];
      print(ifprint2, "col: %i,  obj: %i,  val: %lf,  useful: %lf\n",
    	  object->columnNumber(), infInd_[i],
    	  object->infeasibility(&branchInfo, way), infUseful_[i]);
      print(ifprint2, "    downTot: %lf, downNum: %i, upTot: %lf, upNum: %i\n",
    	  downTotalChange[ind], downNumber[ind],
    	  upTotalChange[ind], upNumber[ind]);
    }
    print(ifprint2, "Before SB Feas\n");
    for (int i = 0; i < feasNum_; ++i) {
      const int ind = feasInd_[i];
      const OsiObject* object = objects[ind];
      print(ifprint2, "colInd: %i,  objInd: %i,  val: %lf,  useful: %lf\n",
    	  object->columnNumber(), feasInd_[i],
    	  object->infeasibility(&branchInfo, way), feasUseful_[i]);
      print(ifprint2, "    downTot: %lf, downNum: %i, upTot: %lf, upNum: %i\n",
    	  downTotalChange[ind], downNumber[ind],
    	  upTotalChange[ind], upNumber[ind]);
    }
  }
#endif  
  
  return infNum_;
}

//-----------------------------------------------------------------------------

void
BM_lp::clear_SB_results()
{
  for (int i = 0; i < objNum_; ++i) {
    sbResult_[i].branchEval = 0;
  }
  bestSbResult_ = NULL;
}

//-----------------------------------------------------------------------------

void
BM_lp::send_one_SB_data(int fixed_size, int changeType, int objInd, int colInd,
			double solval, double bd, int pid)
{
  bm_buf.set_size(fixed_size);
  bm_buf.pack(changeType);
  bm_buf.pack(objInd);
  bm_buf.pack(colInd);
  bm_buf.pack(solval);
  bm_buf.pack(bd);
  send_message(pid, bm_buf);
  BM_SB_result& sbres = sbResult_[objInd];
  sbres.objInd = objInd;
  if (changeType == BM_Var_DownBranch) {
    sbres.varChange[0] = solval - bd;
  } else {
    sbres.varChange[1] = bd - solval;
  }
}

//-----------------------------------------------------------------------------

int
BM_lp::send_data_for_distributed_SB(OsiBranchingInformation& branchInfo,
				    OsiSolverInterface* solver,
				    const int* pids, const int pidNum)
{
  const double * clb = solver->getColLower();
  const double * cub = solver->getColUpper();
  const int numCols = solver->getNumCols();
  bm_buf.clear();
  int tag = BM_StrongBranchRequest;
  bm_buf.pack(tag);
  bm_buf.pack(clb, numCols);
  bm_buf.pack(cub, numCols);
  bm_buf.pack(branchInfo.objectiveValue_);
  bm_buf.pack(branchInfo.cutoff_);
  const int fixed_size = bm_buf.size();
  // We got a few procs to work with, so do distributed strong branching
  int pidCnt = 0;
#ifdef BM_PRINT_DATA
  print(ifprint2, "pidNum: %i,  infNum_: %i,  feasNum_: %i\n",
	pidNum, infNum_, feasNum_);
#endif
  /* For the 0-th object just send out the up-branch, the down branch will be
     processed locally. */
  {
    const int ind = solver->object(infInd_[0])->columnNumber();
    const double val = branchInfo.solution_[ind];
    send_one_SB_data(fixed_size, BM_Var_UpBranch, infInd_[0], ind,
		     val, ceil(val), pids[pidCnt++]);
  }
  int i;
  for (i = 1; i < infNum_ && pidCnt < pidNum; ++i) {
    /* FIXME: think about SOS */
    const int ind = solver->object(infInd_[i])->columnNumber();
    const double val = branchInfo.solution_[ind];
    send_one_SB_data(fixed_size, BM_Var_DownBranch, infInd_[i], ind,
		     val, floor(val), pids[pidCnt++]);
    assert(pidCnt < pidNum);
    send_one_SB_data(fixed_size, BM_Var_UpBranch, infInd_[i], ind,
		     val, ceil(val), pids[pidCnt++]);
  }
  i = 0;
  while (pidCnt < pidNum) {
    print(ifprint2, "running through feas: i: %i   feasInd_[i]: %i\n",
	  i, feasInd_[i]);
    const int ind = solver->object(feasInd_[i])->columnNumber();
    const double lb = branchInfo.lower_[ind];
    const double ub = branchInfo.upper_[ind];
    const double val = branchInfo.solution_[ind];
    if (floor(val+0.5) > lb) { // not at its lb
      send_one_SB_data(fixed_size, BM_Var_DownBranch, feasInd_[i], ind,
		       val, floor(val - 0.5), pids[pidCnt++]);
      if (pidCnt == pidNum) {
	break;
      }
    }
    if (ceil(val-0.5) < ub) { // not at its ub
      send_one_SB_data(fixed_size, BM_Var_UpBranch, feasInd_[i], ind,
		       val, ceil(val + 0.5), pids[pidCnt++]);
    }
    ++i;
  }
  bm_buf.clear();
  return pidCnt;
}

//-----------------------------------------------------------------------------

/* FIXME: this assumes that the solver is the NLP solver. Maybe we should use
   the nlp solver in BM_lp */
void
BM_lp::solve_first_candidate(OsiBranchingInformation& branchInfo,
			     OsiSolverInterface* solver, int downUp)
{
  BM_SB_result& sbres = sbResult_[infInd_[0]];
  sbres.objInd = infInd_[0];
  const int ind = solver->object(infInd_[0])->columnNumber();
  const double val = branchInfo.solution_[ind];
  if (downUp == 0) {
    const double old_bd = solver->getColUpper()[ind];
    solver->setColUpper(ind, floor(val));
    solver->resolve();
    sbres.branchEval |= 1;
    sbres.status[0] =
      (solver->isAbandoned()           ? BCP_Abandoned : 0) |
      (solver->isProvenOptimal()       ? BCP_ProvenOptimal : 0) |
      (solver->isProvenPrimalInfeasible() ? BCP_ProvenPrimalInf : 0);
    sbres.objval[0] =
      (sbres.status[0] & BCP_ProvenOptimal) != 0 ? solver->getObjValue() : 0.0;
    sbres.iter[0] = solver->getIterationCount();
    sbres.varChange[0] = val - floor(val);
    solver->setColUpper(ind, old_bd);
  } else {
    const double old_bd = solver->getColLower()[ind];
    solver->setColLower(ind, ceil(val));
    solver->resolve();
    sbres.branchEval |= 2;
    sbres.status[1] =
      (solver->isAbandoned()           ? BCP_Abandoned : 0) |
      (solver->isProvenOptimal()       ? BCP_ProvenOptimal : 0) |
      (solver->isProvenPrimalInfeasible() ? BCP_ProvenPrimalInf : 0);
    sbres.objval[1] =
      (sbres.status[1] & BCP_ProvenOptimal) != 0 ? solver->getObjValue() : 0.0;
    sbres.iter[1] = solver->getIterationCount();
    sbres.varChange[1] = ceil(val) - val;
    solver->setColLower(ind, old_bd);
  }
}

//-----------------------------------------------------------------------------

void
BM_lp::receive_distributed_SB_result()
{
  int changeType;
  int objInd;
  int colInd;
  int tag;
  bm_buf.clear();
  receive_message(BCP_AnyProcess, bm_buf, BCP_Msg_User);
  bm_buf.unpack(tag);
  assert(tag == BM_StrongBranchResult);
  bm_buf.unpack(changeType);
  bm_buf.unpack(objInd);
  BM_SB_result& sbres = sbResult_[objInd];
  const int field = changeType == BM_Var_UpBranch ? 1 : 0;
  bm_buf.unpack(colInd);
  assert(colInd == sbres.colInd);
  sbres.branchEval |= (changeType == BM_Var_UpBranch ? 2 : 1);
  bm_buf.unpack(sbres.status[field]);
  bm_buf.unpack(sbres.iter[field]);
  bm_buf.unpack(sbres.objval[field]);
}

//-----------------------------------------------------------------------------

bool
BM_lp::isBranchFathomable(int status, double obj)
{
  return ( (status & BCP_ProvenPrimalInf) ||
	   ((status & BCP_ProvenOptimal) && over_ub(obj)) );
}

//-----------------------------------------------------------------------------

int
BM_lp::process_SB_results(OsiBranchingInformation& branchInfo,
			  OsiSolverInterface* solver,
			  Bonmin::BonChooseVariable* choose,
			  OsiBranchingObject*& branchObject)
{
  int i;
#ifdef BM_PRINT_DATA
  print(ifprint2, "SB results Infeas\n");
  for (i = 0; i < infNum_; ++i) {
    const BM_SB_result& sbres = sbResult_[infInd_[i]];
    if (sbres.branchEval == 0) {
      continue;
    }
    print(ifprint2, "eval: %i  col: %i  stati: %i %i,  obj: %f %f\n",
	   sbres.branchEval, sbres.colInd, sbres.status[0], sbres.status[1],
	   sbres.objval[0], sbres.objval[1]);
  }
  print(ifprint2, "SB results Feas\n");
  for (i = 0; i < feasNum_; ++i) {
    const BM_SB_result& sbres = sbResult_[feasInd_[i]];
    if (sbres.branchEval == 0) {
      continue;
    }
    print(ifprint2, "eval: %i  col: %i  stati: %i %i,  obj: %f %f\n",
	   sbres.branchEval, sbres.colInd, sbres.status[0], sbres.status[1],
	   sbres.objval[0], sbres.objval[1]);
  }
#endif
  // First check if we can fathom the node
  int listLen=0; // we want this for the bm_stats
  for (i = 0; i < infNum_; ++i) {
    const BM_SB_result& sbres = sbResult_[infInd_[i]];
    if (sbres.branchEval == 0) {
      continue;
    }
    assert(sbres.branchEval == 3);
    if (isBranchFathomable(sbres.status[0], sbres.objval[0]) &&
	isBranchFathomable(sbres.status[1], sbres.objval[1])) {
      return -2;
    }
    ++listLen;
  }

  // Nope. Let's see if we can fix anything. 
  // 0: nothing fixed   bit 0: fixed but stayd feas  bit 1: fixed and lost feas
  int fixedStat = 0; 
  for (i = 0; i < infNum_; ++i) {
    const BM_SB_result& sbres = sbResult_[infInd_[i]];
    if (sbres.branchEval == 0) {
      continue;
    }
    if (isBranchFathomable(sbres.status[0], sbres.objval[0])) {
      const int colInd = sbres.colInd;
      solver->setColLower(colInd, ceil(branchInfo.solution_[colInd]));
      fixedStat |= 2;
      bm_stats.incNumberFixed();
    }
    if (isBranchFathomable(sbres.status[1], sbres.objval[1])) {
      const int colInd = sbres.colInd;
      solver->setColUpper(colInd, floor(branchInfo.solution_[colInd]));
      fixedStat |= 2;
      bm_stats.incNumberFixed();
    }
  }
  for (i = 0; i < feasNum_; ++i) {
    const BM_SB_result& sbres = sbResult_[feasInd_[i]];
    if (sbres.branchEval == 0) {
      continue;
    }
    if ( (sbres.branchEval & 1) &&
	 isBranchFathomable(sbres.status[0], sbres.objval[0]) ) {
      // Down branch evaluated and fathomable
      const int colInd = sbres.colInd;
      solver->setColLower(colInd, ceil(branchInfo.solution_[colInd] - 0.5));
      fixedStat |= 1;
      bm_stats.incNumberFixed();
    }
    if ( (sbres.branchEval & 2) &&
	 isBranchFathomable(sbres.status[1], sbres.objval[1]) ) {
      // Up branch evaluated and fathomable
      const int colInd = sbres.colInd;
      solver->setColUpper(colInd, floor(branchInfo.solution_[colInd] + 0.5));
      fixedStat |= 1;
      bm_stats.incNumberFixed();
    }
  }
  if ((fixedStat & 2) != 0) {
    return -1;
  }

  // OK, we may have fixed some, but we got to branch, nevertheless.
  // Try to choose something where both sides have finished properly.
  // Note: at this point all stati of the inf pieces is either abandoned,
  // optimal (if primal infeas then we fixed and lost feasibility)
  const bool preferHigh = par.entry(BM_par::PreferHighCombinationInBranching);
  int best = -1;
  int bestWhichWay = 1;
  double bestTrusted = preferHigh ? -COIN_DBL_MAX : COIN_DBL_MAX;
  const double MAXMIN = choose->maxminCrit(&branchInfo);

  for (i = 0; i < infNum_; ++i) {
    const int objInd = infInd_[i];
    const BM_SB_result& sbres = sbResult_[objInd];
    if (sbres.branchEval == 0) {
      // FIXME: Check the pseudocost
      continue;
    }
    if ((sbres.status[0]&BCP_Abandoned) || (sbres.status[1]&BCP_Abandoned)){
      continue;
    }
    double upEst = sbres.objval[1] - branchInfo.objectiveValue_;
    double downEst = sbres.objval[0] - branchInfo.objectiveValue_;
    double value =
      MAXMIN*CoinMin(upEst,downEst) + (1.0-MAXMIN)*CoinMax(upEst,downEst);
    const bool better = ( (  preferHigh && (value > bestTrusted)) ||
			  ( !preferHigh && (value < bestTrusted)) );
    if (better) {
      bestTrusted = value;
      best = i;
      bestWhichWay = upEst > downEst ? 0 : 1;
      // override if there is a preferred way
      const OsiObject* object = solver->object(objInd);
      if (object->preferredWay() >= 0 && object->infeasibility()) {
	bestWhichWay = object->preferredWay();
      }
    }
  }
  if (best == -1) {
    // OMG... In *ALL* evaluated candidates at least one side was abandoned...
    // Oh well, pick the first one not evaluated, or if no such thing then the
    // first where only one side failed, or if no such thing then anything...
    // FIXME: this loop will not be needed when pseudocosts are used
    for (i = 0; i < infNum_; ++i) {
      const BM_SB_result& sbres = sbResult_[infInd_[i]];
      if (sbres.branchEval == 0) {
	best = i;
	break;
      }
    }
    if (best == -1) {
      for (i = 0; i < infNum_; ++i) {
	const BM_SB_result& sbres = sbResult_[infInd_[i]];
	if ((sbres.status[0] & BCP_Abandoned) == 0 ||
	    (sbres.status[1] & BCP_Abandoned) == 0) {
	  // prefer something where at least one side was not abandoned...
	  best = i;
	  break;
	}
	best = i;
      }
    }
  }

  bm_stats.updateStrongBrachingInfo(best, listLen);
  
  // At this point best is not -1, create the branching object
  const OsiObject * object = solver->object(infInd_[best]);
  branchObject = object->createBranch(solver, &branchInfo, bestWhichWay);
  const int ind = object->columnNumber();
  bestSbResult_ = sbResult_ + infInd_[best];
  print(ifprint3,
	"LP: Branch var: %i, val: %lf, obj0: %lf, obj1: %lf, way: %i\n",
	ind, branchInfo.solution_[ind], bestSbResult_->objval[0],
	bestSbResult_->objval[1], bestWhichWay);

  return (fixedStat & 1) != 0 ? -1 : 0;
}

//-----------------------------------------------------------------------------

int
BM_lp::try_to_branch(OsiBranchingInformation& branchInfo,
		     OsiSolverInterface* solver,
		     Bonmin::BonChooseVariable* choose,
		     OsiBranchingObject*& branchObject,
		     bool allowVarFix)
{
  int returnStatus = 0;

  int branchNum; 
  sort_objects(branchInfo, choose, branchNum);

  if (infNum_ == 0) {
    return 0;
  }

  bm_buf.clear();
  bm_buf.pack(branchNum-1);
  send_message(parent(), bm_buf, BCP_Msg_RequestProcessList);
  bm_buf.clear();
  receive_message(parent(), bm_buf, BCP_Msg_ProcessList);
  int* pids = NULL;
  int pidNum;
  bm_buf.unpack(pids, pidNum);

  clear_SB_results();
  // FIXME: This test will have to be changed to >= 1
  if (pidNum >= 0) {
    if (pidNum > 0) {
      pidNum = send_data_for_distributed_SB(branchInfo, solver, pids, pidNum);
    }
    // While the others are working, initialize the result array
    solve_first_candidate(branchInfo, solver, 0);
    if (pidNum == 0) {
      solve_first_candidate(branchInfo, solver, 1);
    }
    while (pidNum > 0) {
      receive_distributed_SB_result();
      --pidNum;
    }
    returnStatus = process_SB_results(branchInfo, solver, choose, branchObject);
    send_pseudo_cost_update(branchInfo);
    
  } else { /* Do something locally */

    returnStatus = BCP_lp_user::try_to_branch(branchInfo, solver, choose,
					      branchObject, allowVarFix);
  }

  delete[] pids;

  return returnStatus;
}

//-----------------------------------------------------------------------------

BCP_branching_decision
BM_lp::bbBranch(OsiBranchingInformation& brInfo,
		BCP_vec<BCP_lp_branching_object*>& cands)
{
  OsiSolverInterface* osi = getLpProblemPointer()->lp_solver;
  Bonmin::OsiTMINLPInterface* nlp =
    dynamic_cast<Bonmin::OsiTMINLPInterface*>(osi);
  assert(nlp);

  nlp->getDblParam(OsiPrimalTolerance, brInfo.integerTolerance_);
    
  BCP_branching_decision retCode;
  OsiBranchingObject* brObj = NULL;

//   static int cnt = 0;
//   print(ifprint, "cnt = %i\n", cnt);
//   ++cnt;

  const int numCols = nlp->getNumCols();
  double* clb_old = new double[numCols];
  double* cub_old = new double[numCols];
  CoinDisjointCopyN(nlp->getColLower(), numCols, clb_old);
  CoinDisjointCopyN(nlp->getColUpper(), numCols, cub_old);

  Ipopt::SmartPtr<Ipopt::OptionsList> options = bonmin_.options();
  int numSB = 0;
  const bool sbIsSet =
    options->GetIntegerValue("number_strong_branch",numSB,"bonmin.");
  int numSBroot = 0;
  const bool sbRootIsSet =
    options->GetIntegerValue("number_strong_branch_root",numSBroot,"bonmin.");

  if (sbRootIsSet && brInfo.depth_ == 0) {
    bonmin_.branchingMethod()->setNumberStrong(numSBroot);
  } else {
    bonmin_.branchingMethod()->setNumberStrong(numSB);
  }

  Bonmin::BonChooseVariable* choose =
    dynamic_cast<Bonmin::BonChooseVariable*>(bonmin_.branchingMethod());
  const int brResult = BM_lp::try_to_branch(brInfo, nlp, choose, brObj, true);

#if 0
  if (choose->goodSolution()) {
    /* yipee! a feasible solution! Check that it's really */
    const double* sol = choose->goodSolution();
    BM_solution* bmsol = new BM_solution;
    //Just copy the solution stored in solver to sol
    double ptol = integerTolerance_;
    for (int i = 0 ; i < numCols ; i++) {
      if (fabs(sol[i]) > ptol)
	bmsol->add_entry(i, sol[i]); 
    }
    bmsol->setObjective(choose->goodObjectiveValue());
    choose->clearGoodSolution();
    send_feasible_solution(bmsol);
    delete bmsol;
  }
#endif

  switch (brResult) {
  case -2:
    // when doing strong branching a candidate has proved that the
    // problem is infeasible
    retCode = BCP_DoNotBranch_Fathomed;
    break;
  case -1:
    // OsiChooseVariable::chooseVariable() returned 2, 3, or 4
    if (!brObj) {
      // just go back and resolve
      retCode = BCP_DoNotBranch;
    } else {
      // otherwise might as well branch. The forced variable is
      // unlikely to jump up one more (though who knows...)
      retCode = BCP_DoBranch;
    }
    break;
  case 0:
    if (!brObj) {
      // nothing got fixed, yet couldn't find something to branch on
      throw BCP_fatal_error("BM: Couldn't branch!\n");
    }
    // we've got a branching object
    retCode = BCP_DoBranch;
    break;
  default:
    throw BCP_fatal_error("\
BM: BCP_lp_user::try_to_branch returned with unknown return code.\n");
  }

  if (brResult < 0) {
    // If there were some fixings (brResult < 0) then propagate them
    // where needed

    // FIXME: This is not nice. Meddling w/ BCP internal data. The BCP
    // user interface should provide a way to change bounds regardless
    // whether branching is asked for or not.
    const double* clb = nlp->getColLower();
    const double* cub = nlp->getColUpper();

    BCP_vec<BCP_var*>& vars = getLpProblemPointer()->node->vars;
    for (int i = 0; i < numCols; ++i) {
      if (clb_old[i] != clb[i] || cub_old[i] != cub[i]) {
	vars[i]->change_bounds(clb[i], cub[i]);
	nlp->setColBounds(i, clb[i], cub[i]);
      }
    }
  }

  if (retCode == BCP_DoBranch) {
    // all possibilities are 2-way branches
    int order[2] = {0, 1};
    // Now interpret the result (at this point we must have a brObj
    OsiIntegerBranchingObject* intBrObj =
      dynamic_cast<OsiIntegerBranchingObject*>(brObj);
    if (intBrObj) {
      if (intBrObj->firstBranch() == 1) {
	order[0] = 1;
	order[1] = 0;
      }
      BCP_lp_integer_branching_object o(intBrObj);
      cands.push_back(new BCP_lp_branching_object(o, order));
      if (bestSbResult_) {
	BCP_vec<double> lb(2, 0.0);
	lb[0] = bestSbResult_->objval[order[0]];
	lb[1] = bestSbResult_->objval[order[1]];
	BCP_vec<int> tc(2, 0);
	tc[0] = bestSbResult_->status[order[0]];
	tc[1] = bestSbResult_->status[order[1]];
	cands.back()->set_presolve_result(lb, tc);
      }
      if (par.entry(BM_par::PrintBranchingInfo)) {
	print(ifprint2, "BM_lp: branching on variable %i   value: %f\n",
	       intBrObj->originalObject()->columnNumber(),
	       intBrObj->value());
      }
    }
    OsiSOSBranchingObject* sosBrObj =
      dynamic_cast<OsiSOSBranchingObject*>(brObj);
    if (sosBrObj) {
      if (sosBrObj->firstBranch() == 1) {
	order[0] = 1;
	order[1] = 0;
      }
      BCP_lp_sos_branching_object o(sosBrObj);
      cands.push_back(new BCP_lp_branching_object(nlp, o, order));
      if (bestSbResult_) {
	BCP_vec<double> lb(2, 0.0);
	lb[0] = bestSbResult_->objval[order[0]];
	lb[1] = bestSbResult_->objval[order[1]];
	BCP_vec<int> tc(2, 0);
	tc[0] = bestSbResult_->status[order[0]];
	tc[1] = bestSbResult_->status[order[1]];
	cands.back()->set_presolve_result(lb, tc);
      }
      if (par.entry(BM_par::PrintBranchingInfo)) {
	print(ifprint2, "BM_lp: branching on SOS %i   value: %f\n",
	       sosBrObj->originalObject()->columnNumber(),
	       sosBrObj->value());
      }
    }
  }

  delete brObj;
  delete[] clb_old;
  delete[] cub_old;
  return retCode;
}

/****************************************************************************/

void
BM_lp::set_user_data_for_children(BCP_presolved_lp_brobj* best, 
                                  const int selected)
{
    BM_node* data = NULL;
    data = new BM_node;
    data->numNlpFailed_ = numNlpFailed_;
    best->user_data()[0] = data;
    data = new BM_node;
    data->numNlpFailed_ = numNlpFailed_;
    best->user_data()[1] = data;
}

//#############################################################################

void
BM_lp::process_message(BCP_buffer& buf)
{
  int msgtag;
  buf.unpack(msgtag);
  assert(msgtag == BM_StrongBranchRequest);

  Bonmin::OsiTMINLPInterface& nlp = *bonmin_.nonlinearSolver();

  int numCols;
  double* clb;
  double* cub;
  double objvalOrig;
  double cutoff;
  buf.unpack(clb, numCols);
  assert(numCols == nlp.getNumCols());
  buf.unpack(cub, numCols);
  assert(numCols == nlp.getNumCols());
  buf.unpack(objvalOrig);
  buf.unpack(cutoff);
  int changeType;
  int objInd;
  int colInd;
  double solval;
  double bd;
  buf.unpack(changeType);
  buf.unpack(objInd);
  buf.unpack(colInd);
  buf.unpack(solval);
  buf.unpack(bd);
  nlp.setColLower(clb);
  nlp.setColUpper(cub);
  switch (changeType) {
  case BM_Var_DownBranch:
    nlp.setColUpper(colInd, bd);
    break;
  case BM_Var_UpBranch:
    nlp.setColLower(colInd, bd);
    break;
  }
  nlp.resolve();
  int status = 
    (nlp.isAbandoned()              ? BCP_Abandoned : 0) |
    (nlp.isProvenOptimal()          ? BCP_ProvenOptimal : 0) |
    (nlp.isProvenPrimalInfeasible() ? BCP_ProvenPrimalInf : 0);
  double objval = (status & BCP_ProvenOptimal) != 0 ? nlp.getObjValue() : 0.0;
  int iter = nlp.getIterationCount();

  bm_buf.clear();
  msgtag = BM_StrongBranchResult;
  bm_buf.pack(msgtag);
  bm_buf.pack(changeType);
  bm_buf.pack(objInd);
  bm_buf.pack(colInd);
  bm_buf.pack(status);
  bm_buf.pack(iter);
  bm_buf.pack(objval);
  send_message(buf.sender(), bm_buf, BCP_Msg_User);

  bm_buf.clear();
  bm_buf.pack(changeType);
  bm_buf.pack(objvalOrig);
  bm_buf.pack(cutoff);
  bm_buf.pack(objInd);
  bm_buf.pack(colInd);
  bm_buf.pack(solval);
  bm_buf.pack(bd);
  bm_buf.pack(status);
  bm_buf.pack(iter);
  bm_buf.pack(objval);
  send_message(parent(), bm_buf, BCP_Msg_SBnodeFinished);

  delete[] clb;
  delete[] cub;

  bm_stats.incNumberSbSolves();
}

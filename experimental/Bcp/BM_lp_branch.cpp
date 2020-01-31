// (C) Copyright International Business Machines Corporation 2006, 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Laszlo Ladanyi, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University

#include "CoinHelperFunctions.hpp"
#include "BCP_lp_node.hpp"
#include "BCP_lp.hpp"
#include "BCP_lp_functions.hpp"
#include "BM.hpp"
#include "BonChooseVariable.hpp"
#include "BonCurvBranchingSolver.hpp"
#include "BonQpBranchingSolver.hpp"
#include "BonLpBranchingSolver.hpp"
#include "BonOsiTMINLPInterface.hpp"
#include "BonIpoptWarmStart.hpp"

#ifndef BM_DEBUG_PRINT
#define BM_DEBUG_PRINT 0
#endif

static bool ifprint = true;
static bool ifprint2 = false;

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

  const bool disregardPriorities = par.entry(BM_par::DisregardPriorities);
  const bool usePseudoCosts = par.entry(BM_par::UsePseudoCosts);

  for (int i = 0; i < objNum_; ++i) {
    const int ind = objInd_[i];
    const OsiObject* object = objects[ind];
    double value = object->infeasibility(&branchInfo, way);
    if (value  > 0.0) {
      if (value >= 1e50) { // infeasible
	return -1;
      }
      if (! disregardPriorities) {
	int priorityLevel = object->priority();
	if (lastPriority < priorityLevel) {
	  // sort the entries based on their usefulness
	  if (infBlockStart < infNum_) {
	    if (par.entry(BM_par::DecreasingSortInSetupList)) {
	      CoinSort_2(infUseful_ + infBlockStart, infUseful_ + infNum_,
			 infInd_ + infBlockStart,
			 CoinFirstGreater_2<double,int>());
	    } else {
	      CoinSort_2(infUseful_ + infBlockStart, infUseful_ + infNum_,
			 infInd_ + infBlockStart);
	    }
	    infBlockStart = infNum_;
	  }
	  lastPriority = priorityLevel;
	}
      }
      double dummy;
      infInd_[infNum_] = ind;
      if (usePseudoCosts) {
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
      if (! disregardPriorities) {
	int priorityLevel = object->priority();
	if (lastPriority < priorityLevel) {
	  // sort the entries based on their usefulness
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
	  lastPriority = priorityLevel;
	}
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

#if (BM_DEBUG_PRINT != 0)
  const double t = CoinWallclockTime();
  printf("LP %.3f: node: %i  depth: %i  obj: %f  infNum: %i  feasNum: %i  soltime: %.3f\n",
	 t-start_time(), current_index(), current_level(),
	 branchInfo.objectiveValue_, 
	 infNum_, feasNum_, t-node_start_time);
  node_start_time = t;
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
BM_lp::collect_branch_data(OsiBranchingInformation& branchInfo,
			   OsiSolverInterface* solver,
			   const int branchNum,
			   BM_BranchData* branchData)
{
  int i;
  int b = 0;
  for (i = 0; i < infNum_; ++i) {
    /* FIXME: think about SOS */
    const int objInd = infInd_[i];
    const int colInd = solver->object(objInd)->columnNumber();
    const double val = branchInfo.solution_[colInd];
    branchData[b].changeType = BM_Var_DownBranch;
    branchData[b].objInd = objInd;
    branchData[b].colInd = colInd;
    branchData[b].solval = val;
    branchData[b].bd = floor(val);
    BM_SB_result& sbres = sbResult_[objInd];
    sbres.objInd = objInd;
    sbres.varChange[0] = val - floor(val);
    ++b;
    if (b == branchNum) {
      return;
    }
    branchData[b].changeType = BM_Var_UpBranch;
    branchData[b].objInd = objInd;
    branchData[b].colInd = colInd;
    branchData[b].solval = val;
    branchData[b].bd = ceil(val);
    sbres.varChange[1] = ceil(val) - val;
    ++b;
    if (b == branchNum) {
      return;
    }
  }
  for (i = 0; i < feasNum_; ++i) {
    const int objInd = feasInd_[i];
    const int colInd = solver->object(objInd)->columnNumber();
    const double val = branchInfo.solution_[colInd];
    const double lb = branchInfo.lower_[colInd];
    const double ub = branchInfo.upper_[colInd];
    if (floor(val+0.5) > lb) { // not at its lb
      branchData[b].changeType = BM_Var_DownBranch;
      branchData[b].objInd = objInd;
      branchData[b].colInd = colInd;
      branchData[b].solval = val;
      branchData[b].bd = floor(val - 0.5);
      ++b;
      if (b == branchNum) {
	return;
      }
    }
    if (ceil(val-0.5) < ub) { // not at its ub
      branchData[b].changeType = BM_Var_UpBranch;
      branchData[b].objInd = objInd;
      branchData[b].colInd = colInd;
      branchData[b].solval = val;
      branchData[b].bd = ceil(val + 0.5);
      ++b;
      if (b == branchNum) {
	return;
      }
    }
  }
}

//-----------------------------------------------------------------------------

void
BM_solve_branches(OsiSolverInterface* solver, const CoinWarmStart* cws,
		  const int numBranch, BM_BranchData* bD)
{
  for (int i = 0; i < numBranch; ++i) {
    double t = CoinWallclockTime();
    const int ind = bD[i].colInd;
    const int field = bD[i].changeType == BM_Var_UpBranch ? 1 : 0;
    const double old_lb = solver->getColLower()[ind];
    const double old_ub = solver->getColUpper()[ind];
    if (field == 0) {
      solver->setColUpper(ind, bD[i].bd);
    } else {
      solver->setColLower(ind, bD[i].bd);
    }
    if (cws) {
      solver->setWarmStart(cws);
      solver->resolve();
    } else {
      solver->initialSolve();
    }
    bD[i].status =
      (solver->isAbandoned()              ? BCP_Abandoned : 0) |
      (solver->isProvenOptimal()          ? BCP_ProvenOptimal : 0) |
      (solver->isProvenPrimalInfeasible() ? BCP_ProvenPrimalInf : 0);
    bD[i].objval =
      (bD[i].status & BCP_ProvenOptimal) != 0 ? solver->getObjValue() : 0.0;
    bD[i].iter = solver->getIterationCount();
    solver->setColBounds(ind, old_lb, old_ub);
    bD[i].time = CoinWallclockTime() - t;
  }
}

//-----------------------------------------------------------------------------

void
BM_register_branch_results(const int numBranch, const BM_BranchData* bD,
			   BM_SB_result* sbResults)
{
  for (int i = 0; i < numBranch; ++i) {
    const int field = bD[i].changeType == BM_Var_UpBranch ? 1 : 0;
    BM_SB_result& sbres = sbResults[bD[i].objInd];
    sbres.objInd           = bD[i].objInd;
    sbres.branchEval      |= field == 0 ? 1 : 2;
    sbres.status[field]    = bD[i].status;
    sbres.objval[field]    = bD[i].objval;
    sbres.iter[field]      = bD[i].iter;
    sbres.varChange[field] = fabs(bD[i].solval - bD[i].bd);
  }
}

//-----------------------------------------------------------------------------

void
BM_lp::do_distributed_SB(OsiBranchingInformation& branchInfo,
			 OsiSolverInterface* solver,
			 const CoinWarmStart* cws,
			 const int branchNum,
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
  bool has_ws = cws != NULL;
  bm_buf.pack(has_ws);
  if (has_ws) {
    BCP_lp_prob* p = getLpProblemPointer();
    CoinWarmStart* cws_tmp = cws->clone(); // the next call destroys it
    BCP_warmstart* ws = cws ? BCP_lp_convert_CoinWarmStart(*p, cws_tmp) : NULL;
    BCP_pack_warmstart(ws, bm_buf);
    delete ws;
  }
  
  const int fixed_size = bm_buf.size();

  // collect what we'll need to send off data
  BM_BranchData* branchData = new BM_BranchData[branchNum];
  collect_branch_data(branchInfo, solver, branchNum, branchData);

  // We have branchNum branches to process on pidNum+1 (the last is the local
  // process) processes.
  int branchLeft = branchNum;
  BM_BranchData* bD = branchData;
  for (int pidLeft = pidNum; pidLeft > 0; --pidLeft) {
    int numSend = branchLeft / pidLeft;
    if (numSend * pidLeft < branchLeft) {
      ++numSend;
    }
    bm_buf.set_size(fixed_size);
    // Now pack where we are in branchData and pack numSend branches
    int location = bD - branchData;
    bm_buf.pack(location);
    bm_buf.pack(numSend);
    for (int s = 0; s < numSend; ++s) {
      bm_buf.pack(bD[s].changeType);
      bm_buf.pack(bD[s].objInd);
      bm_buf.pack(bD[s].colInd);
      bm_buf.pack(bD[s].solval);
      bm_buf.pack(bD[s].bd);
    }
    send_message(pids[pidLeft-1], bm_buf);
    bD += numSend;
    branchLeft -= numSend;
  }
  assert(branchNum/(pidNum+1) == branchLeft);

  // Process the leftover branches locally
  /* FIXME: this assumes that the solver is the NLP solver. Maybe we should use
     the nlp solver in BM_lp */
  BM_solve_branches(solver, cws, branchLeft, bD);
  bm_stats.incNumberSbSolves(branchLeft);
  BM_register_branch_results(branchLeft, bD, sbResult_);

  // Receive the results from the other processes
  int numResults = branchLeft;
  while (numResults < branchNum) {
    bm_buf.clear();
    receive_message(BCP_AnyProcess, bm_buf, BCP_Msg_User);
    bm_buf.unpack(tag);
    assert(tag == BM_StrongBranchResult);
    int location;
    int numRes;
    bm_buf.unpack(location);
    bm_buf.unpack(numRes);
    BM_BranchData* bD = branchData + location;
    for (int i = 0; i < numRes; ++i) {
      bm_buf.unpack(bD[i].status);
      bm_buf.unpack(bD[i].objval);
      bm_buf.unpack(bD[i].iter);
      bm_buf.unpack(bD[i].time);
    }
    BM_register_branch_results(numRes, bD, sbResult_);
    numResults += numRes;
  }
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
  int evaluated = 0;
  int i;
#if defined(DEBUG_PRINT)
  const double t = CoinWallclockTime()-start_time();
#endif

#if 0 // defined(DEBUG_PRINT)
  const double t = CoinWallclockTime()-start_time();
  for (i = 0; i < infNum_; ++i) {
    const BM_SB_result& sbres = sbResult_[infInd_[i]];
    if (sbres.branchEval == 0) {
      continue;
    }
    printf("LP %.3f: SB: node: %i  inf col: %i  stati: %i %i,  obj: %f %f  time: %.3f %.3f\n",
	   t, current_index(), sbres.colInd, 
	   sbres.status[0], sbres.status[1],
	   sbres.objval[0], sbres.objval[1], sbres.time[0], sbres.time[1]);
  }
  for (i = 0; i < feasNum_; ++i) {
    const BM_SB_result& sbres = sbResult_[feasInd_[i]];
    if (sbres.branchEval == 0) {
      continue;
    }
    printf("LP %.3f: SB: node: %i  feas col: %i  stati: %i %i,  obj: %f %f  time: %.3f %.3f\n",
	   t, current_index(), sbres.colInd, 
	   sbres.status[0], sbres.status[1],
	   sbres.objval[0], sbres.objval[1], sbres.time[0], sbres.time[1]);
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
#if (BM_DEBUG_PRINT != 0)
      const double wallclock = CoinWallclockTime();
#if 0
      printf("LP %.3f: SBres: node: %i  FATHOM  inf/eval/cand: time: %.3f\n",
	     wallclock-start_time(), current_index(),
	     infNum_, 
	     wallclock-node_start_time);
#else 
      printf("LP %.3f: SBres: node: %i  FATHOM  time: %.3f\n",
	     wallclock-start_time(), current_index(),
	     wallclock-node_start_time);
#endif
#endif
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
    ++evaluated;
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
    ++evaluated;
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
#if (BM_DEBUG_PRINT != 0)
    printf("LP %.3f: SBres: node: %i  RESOLVE  time: %.3f\n",
	   t - start_time(), current_index(), t-node_start_time);
    node_start_time = t;
#endif
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
  bestSbResult_ = sbResult_ + infInd_[best];
#if (BM_DEBUG_PRINT != 0)
  const int ind = object->columnNumber();
  printf("LP %.3f: SBres: node: %i  depth: %i  BRANCH  time: %.3f  evaluated: %i  bvar: %i  val: %f  obj0: %f  obj1: %f  way: %i\n",
	 t-start_time(), current_index(), current_level(), t-node_start_time, 
	 evaluated, ind, branchInfo.solution_[ind], 
	 bestSbResult_->objval[0], bestSbResult_->objval[1], bestWhichWay);
  node_start_time = t;
#endif

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

  const bool do_distributed_branching = true;

  if (do_distributed_branching) {

    bm_buf.clear();
    bm_buf.pack(branchNum-1);
    send_message(parent(), bm_buf, BCP_Msg_RequestProcessList);
    bm_buf.clear();
    receive_message(parent(), bm_buf, BCP_Msg_ProcessList);
    int* pids = NULL;
    int pidNum;
    bm_buf.unpack(pids, pidNum);

    clear_SB_results();

    // Get the warmstart information
    CoinWarmStart* cws = solver->getWarmStart();
    Bonmin::IpoptWarmStart* iws = dynamic_cast<Bonmin::IpoptWarmStart*>(cws);
    if (iws && iws->empty()) {
      delete cws;
      cws = NULL;
    }

    const int nodeIndex = current_index();
    if (nodeIndex == 0) {
      // we are in the root. We want to process all possible branches upto the
      // parameter value.
      branchNum = CoinMin(branchNum, par.entry(BM_par::SBNumBranchesInRoot));
    } else {
      if (nodeIndex < par.entry(BM_par::SBMaxLevel)) {
	branchNum = CoinMin(branchNum,
			    CoinMax(pidNum + 1,
				    par.entry(BM_par::SBNumBranchesInTree)));
      } else {
	branchNum = pidNum > 0 ? pidNum + 1 : 0;
      }
    }

    if (branchNum > 0) {
      do_distributed_SB(branchInfo, solver, cws, branchNum, pids, pidNum);
      returnStatus = process_SB_results(branchInfo, solver, choose,
					branchObject);
      send_pseudo_cost_update(branchInfo);
    } else {
      // We are too deep in the tree, just do something locally (like
      // pseudocosts...)
      returnStatus = BCP_lp_user::try_to_branch(branchInfo, solver, choose,
						branchObject, allowVarFix);
    }

    delete[] pids;
    delete cws;

  } else { /* ! do_distributed_branching  ===> Do something locally */
    
  }

  return returnStatus;
}

//-----------------------------------------------------------------------------

BCP_branching_decision
BM_lp::bbBranch(OsiBranchingInformation& brInfo,
		BCP_vec<BCP_lp_branching_object*>& cands)
{
  BCP_lp_prob* p = getLpProblemPointer();
  OsiSolverInterface* osi = p->lp_solver;
  Bonmin::OsiTMINLPInterface* nlp =
    dynamic_cast<Bonmin::OsiTMINLPInterface*>(osi);
  assert(nlp);

  nlp->getDblParam(OsiPrimalTolerance, brInfo.integerTolerance_);
    
  BCP_branching_decision retCode;
  OsiBranchingObject* brObj = NULL;

  const int numCols = nlp->getNumCols();
  double* clb_old = new double[numCols];
  double* cub_old = new double[numCols];
  CoinDisjointCopyN(nlp->getColLower(), numCols, clb_old);
  CoinDisjointCopyN(nlp->getColUpper(), numCols, cub_old);

  Ipopt::SmartPtr<Ipopt::OptionsList> options = bonmin_.options();
  int numSB = 0;
  //const bool sbIsSet =
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

  Bonmin::OsiTMINLPInterface* nlp = bonmin_.nonlinearSolver();

  int numCols;
  double* clb;
  double* cub;
  double objvalOrig;
  double cutoff;
  buf.unpack(clb, numCols);
  assert(numCols == nlp->getNumCols());
  buf.unpack(cub, numCols);
  assert(numCols == nlp->getNumCols());
  buf.unpack(objvalOrig);
  buf.unpack(cutoff);
  bool has_ws;
  buf.unpack(has_ws);
  BCP_warmstart* ws = has_ws ? BCP_unpack_warmstart(buf) : NULL;
  CoinWarmStart* cws = ws  ? ws->convert_to_CoinWarmStart() : NULL;
  int location; // just a marker that we'll send back
  int numBranch;
  buf.unpack(location);
  buf.unpack(numBranch);

  int i;
  BM_BranchData* branchData = new BM_BranchData[numBranch];
  for (i = 0; i < numBranch; ++i) {
    buf.unpack(branchData[i].changeType);
    buf.unpack(branchData[i].objInd);
    buf.unpack(branchData[i].colInd);
    buf.unpack(branchData[i].solval);
    buf.unpack(branchData[i].bd);
  }

  nlp->setColLower(clb);
  nlp->setColUpper(cub);
  BM_solve_branches(nlp, cws, numBranch, branchData);

  bm_buf.clear();
  msgtag = BM_StrongBranchResult;
  bm_buf.pack(msgtag);
  bm_buf.pack(location);
  bm_buf.pack(numBranch);
  for (i = 0; i < numBranch; ++i) {
    const BM_BranchData& bD = branchData[i];
    bm_buf.pack(bD.status);
    bm_buf.pack(bD.objval);
    bm_buf.pack(bD.iter);
    bm_buf.pack(bD.time);
  }
  send_message(buf.sender(), bm_buf, BCP_Msg_User);

  bm_buf.clear();
  send_message(parent(), bm_buf, BCP_Msg_SBnodeFinished);

  delete[] branchData;
  delete cws;
  delete ws;
  delete[] clb;
  delete[] cub;

  bm_stats.incNumberSbSolves(numBranch);
}

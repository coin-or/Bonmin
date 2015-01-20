// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  05/26/2005

#include "BonOaNlpOptim.hpp"
#include "OsiAuxInfo.hpp"
#include "CbcModel.hpp"
#include "BonBabInfos.hpp"
#include "BonCbc.hpp"

namespace Bonmin
{
   static const char * txt_id = "NLP relax. for OA";

/// Default constructor
  OaNlpOptim::OaNlpOptim(OsiTMINLPInterface * si,
      int maxDepth, bool addOnlyViolated, bool global)
      :
      CglCutGenerator(),
      nlp_(si),
      maxDepth_(maxDepth),
      nSolve_(0),
      addOnlyViolated_(addOnlyViolated),
      global_(global)
  {
    handler_ = new CoinMessageHandler();
    handler_ -> setLogLevel(1);
    messages_ = OaMessages();
  }

  OaNlpOptim::OaNlpOptim(BabSetupBase &b):
      CglCutGenerator(),
      nlp_(b.nonlinearSolver()),
      maxDepth_(1000),
      nSolve_(0)
  {
    int ivalue;
    b.options()->GetEnumValue("add_only_violated_oa", ivalue,b.prefix());
    addOnlyViolated_ = ivalue;
    b.options()->GetEnumValue("oa_cuts_scope", ivalue,b.prefix());
    global_ = ivalue;

    b.options()->GetIntegerValue("nlp_solve_max_depth", maxDepth_,b.prefix());
    b.options()->GetNumericValue("nlp_solves_per_depth", solves_per_level_,b.prefix());
    handler_ = new CoinMessageHandler();
    handler_ -> setLogLevel(1);
    messages_ = OaMessages();
  }
/// Assign an OsiTMINLPInterface
  void
  OaNlpOptim::assignInterface(OsiTMINLPInterface * si)

  {
    nlp_ = si;
  }
/// cut generation method
  void
  OaNlpOptim::generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
      const CglTreeInfo info) {
    if (nlp_ == NULL) {
      CoinError("Error in cut generator for outer approximation no ipopt NLP assigned", "generateCuts", "OaNlpOptim");
    }

    int numcols = nlp_->getNumCols();

    //Get the continuous solution
    //const double *colsol = si.getColSolution();
    //Check for integer feasibility
   if(!info.inTree || info.pass > 0) return;
#if 1
    BabInfo * babInfo = dynamic_cast<BabInfo *> (si.getAuxiliaryInfo());
    assert(babInfo);
    assert(babInfo->babPtr());
    const CbcNode * node = babInfo->babPtr()->model().currentNode();
    int level = (node == NULL) ? 0 : babInfo->babPtr()->model().currentNode()->depth();
    if (info.level > maxDepth_)
      return;
    if(solves_per_level_ < 1e10){
      double rand = CoinDrand48();
      double score = pow(2.,-level)*solves_per_level_;
      //printf("depth %i, score %g , rand %g\n", level, score, rand);
      if (score <= rand)
        return;
    }
#endif
    //Fix the variable which have to be fixed, after having saved the bounds
    double * saveColLb = new double[numcols];
    double * saveColUb = new double[numcols];
    CoinCopyN(nlp_->getColLower(), numcols , saveColLb);
    CoinCopyN(nlp_->getColUpper(), numcols , saveColUb);
    for (int i = 0; i < numcols ; i++) {
      if (nlp_->isInteger(i)) {
        nlp_->setColBounds(i,si.getColLower()[i],si.getColUpper()[i]);
      }
    }

    //Now solve the NLP get the cuts, reset bounds and get out

    //  nlp_->turnOnIpoptOutput();
    nSolve_++;
    nlp_->resolve(txt_id);
    const double * violatedPoint = (addOnlyViolated_)? si.getColSolution():
        NULL;
    nlp_->getOuterApproximation(cs, 1, violatedPoint,global_);
    if (nlp_->isProvenOptimal()) {
      handler_->message(LP_ERROR,messages_)
      <<nlp_->getObjValue()-si.getObjValue()<<CoinMessageEol;
      bool feasible = 1;
      const double * colsol2 = nlp_->getColSolution();
      for (int i = 0 ; i < numcols && feasible; i++) {
        if (nlp_->isInteger(i)) {
          if (fabs(colsol2[i] - floor(colsol2[i] + 0.5) ) > 1e-07)
            feasible = 0;
        }
      }
      if (feasible ) {
#if 1
        // Also store into solver
        OsiAuxInfo * auxInfo = si.getAuxiliaryInfo();
        OsiBabSolver * auxiliaryInfo = dynamic_cast<OsiBabSolver *> (auxInfo);
        if (auxiliaryInfo) {
          double * lpSolution = new double[numcols + 1];
          CoinCopyN(colsol2, numcols, lpSolution);
          lpSolution[numcols] = nlp_->getObjValue();
          auxiliaryInfo->setSolution(lpSolution, numcols + 1, lpSolution[numcols]);
          delete [] lpSolution;
        }
        else
          fprintf(stderr,"No auxiliary info in nlp solve!\n");
#endif

      }
    }
    else if (nlp_->isAbandoned() || nlp_->isIterationLimitReached()) {
      throw CoinError("Unsolved NLP ... exit", "generateCuts", "OaNlpOptim");
    }
    else {
      //       //Add an infeasibility local constraint
      //       CoinPackedVector v;
      //       double rhs = 1.;
      //       for(int i = 0; i < numcols ; i++)
      // 	{
      // 	  if(nlp_->isInteger(i) && (si.getColUpper()[i] - si.getColLower()[i] < 0.9))
      // 	    {
      // 	      double value = floor(colsol[i] + 0.5);
      // 	      assert(fabs(colsol[i]-value)<1e-8 && value >= -1e-08 && value <= 1+ 1e08);
      // 	      v.insert(i, -(2*value - 1));
      // 	      rhs -= value;
      // 	    }
      // 	}
      //       OsiRowCut cut;
      //       cut.setRow(v);
      //       cut.setLb(rhs);
      //       cut.setUb(1e300);
      //       cut.setGloballyValid();
      //       cs.insert(cut);
    }
    for (int i = 0; i < numcols ; i++) {
      if (nlp_->isInteger(i)) {
        nlp_->setColBounds(i,saveColLb[i],saveColUb[i]);
      }
    }
#if 0
    nlp_->deleteLastRows(numberCuts);
#endif
    delete [] saveColLb;
    delete [] saveColUb;
  }

  void
  OaNlpOptim::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("NLP solves in hybrid algorithm (B-Hyb)", RegisteredOptions::BonminCategory);
    roptions->AddLowerBoundedIntegerOption("nlp_solve_frequency",
        "Specify the frequency (in terms of nodes) at which NLP relaxations are solved in B-Hyb.",
        0,10,
        "A frequency of 0 amounts to to never solve the NLP relaxation.");
    roptions->setOptionExtraInfo("nlp_solve_frequency",1);
    roptions->AddLowerBoundedIntegerOption("nlp_solve_max_depth",
        "Set maximum depth in the tree at which NLP relaxations are solved in B-Hyb.",
        0,10,
        "A depth of 0 amounts to to never solve the NLP relaxation.");
    roptions->setOptionExtraInfo("nlp_solve_max_depth",1);
    roptions->AddLowerBoundedNumberOption("nlp_solves_per_depth",
        "Set average number of nodes in the tree at which NLP relaxations are solved in B-Hyb for each depth.",
        0.,false,1e100);
    roptions->setOptionExtraInfo("nlp_solves_per_depth",1);
  }


}/* End namespace Bonmin. */

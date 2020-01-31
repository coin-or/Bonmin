// (C) Copyright International Business Machines 2006-2011
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date : 12/26/2006

#include "BonOaFeasChecker.hpp"
#include "BonminConfig.h"

#include "OsiAuxInfo.hpp"
#include "BonSolverHelp.hpp"

namespace Bonmin
{
   static const char * txt_id = "check integer sol.";

  /// New usefull constructor
  OaFeasibilityChecker::OaFeasibilityChecker(BabSetupBase &b):
      OaDecompositionBase(b, false, true),
      cut_count_(0)
  {
    int ival;
    b.options()->GetEnumValue("feas_check_cut_types", ival, b.prefix()); 
    type_ = CutsTypes(ival);
    b.options()->GetEnumValue("feas_check_discard_policy", ival, b.prefix()); 
    pol_ = CutsPolicies(ival);
    b.options()->GetIntegerValue("generate_benders_after_so_many_oa", ival, b.prefix());
    maximum_oa_cuts_ = static_cast<unsigned int>(ival);
  }
  OaFeasibilityChecker ::~OaFeasibilityChecker ()
  {}

  /// OaDecomposition method
  double
  OaFeasibilityChecker::performOa(OsiCuts & cs, solverManip &lpManip,
      BabInfo * babInfo, double &cutoff,const CglTreeInfo & info) const
  {
    bool isInteger = true;
    bool feasible = 1;

    OsiSolverInterface * lp = lpManip.si();
    OsiBranchingInformation branch_info(lp,false);
    //int numcols = lp->getNumCols();
    double milpBound = -COIN_DBL_MAX;
    int numberPasses = 0;
    double * nlpSol =  NULL;
    int numberCutsBefore = cs.sizeRowCuts();
   
    while (isInteger && feasible ) {
      numberPasses++;

      //setup the nlp

      //Fix the variable which have to be fixed, after having saved the bounds
      double * colsol = const_cast<double *>(lp->getColSolution());
      branch_info.solution_ = colsol;
      fixIntegers(*nlp_,branch_info, parameters_.cbcIntegerTolerance_,objects_, nObjects_);


      //Now solve the NLP get the cuts, and intall them in the local LP
      nlp_->resolve(txt_id);
      if (post_nlp_solve(babInfo, cutoff)) {
        //nlp solved and feasible
        // Update the cutoff
        double ub = nlp_->getObjValue();
        cutoff = ub > 0 ? ub *(1 - parameters_.cbcCutoffIncrement_) : ub*(1 + parameters_.cbcCutoffIncrement_);
        // Update the lp solver cutoff
        lp->setDblParam(OsiDualObjectiveLimit, cutoff);
      }
      // Get the cuts outer approximation at the current point

      nlpSol = const_cast<double *>(nlp_->getColSolution());

      const double * toCut = (parameter().addOnlyViolated_)?
                             colsol:NULL;
      if(cut_count_ <= maximum_oa_cuts_ && type_ == OA)
        nlp_->getOuterApproximation(cs, nlpSol, 1, toCut,
                                    true);
      else {//if (type_ == Benders)
        nlp_->getBendersCut(cs, parameter().global_);
      }
      if(pol_ == DetectCycles)
        nlp_->getBendersCut(savedCuts_, parameter().global_);

      int numberCuts = cs.sizeRowCuts() - numberCutsBefore;
      cut_count_ += numberCuts;
      if (numberCuts > 0)
        installCuts(*lp, cs, numberCuts);

      lp->resolve();
      double objvalue = lp->getObjValue();
      //milpBound = max(milpBound, lp->getObjValue());
      feasible = (lp->isProvenOptimal() &&
          !lp->isDualObjectiveLimitReached() && (objvalue<cutoff)) ;
      //if value of integers are unchanged then we have to get out
      bool changed = true;//if lp is infeasible we don't have to check anything
      isInteger = 0;
      //	  if(!fixed)//fathom on bounds
      //           milpBound = 1e200;
      if (feasible) {
        changed = isDifferentOnIntegers(*nlp_, objects_, nObjects_,
                                        0.1,
                                        nlp_->getColSolution(), lp->getColSolution());
      }
      if (changed) {
       branch_info.solution_ = lp->getColSolution();
       isInteger = integerFeasible(*lp,branch_info, parameters_.cbcIntegerTolerance_,
                                     objects_, nObjects_);
      }
      else {
        isInteger = 0;
        //	  if(!fixed)//fathom on bounds
         milpBound = 1e200;
      }
#ifdef OA_DEBUG
      printf("Obj value after cuts %g, %d rows\n",lp->getObjValue(),
          numberCuts) ;
#endif
    }
    int num_cuts_now = cs.sizeRowCuts();
    if(pol_ == KeepAll){
      for(int i = numberCutsBefore ; i < num_cuts_now ; i++){
        cs.rowCut(i).setEffectiveness(99.9e99);
      }
    }

#ifdef OA_DEBUG
    debug_.printEndOfProcedureDebugMessage(cs, true, cutoff, milpBound, isInteger, feasible, std::cout);
    std::cout<<"milpBound found: "<<milpBound<<std::endl;
#endif
    return milpBound;
  }


  /** Register OA feasibility checker  options.*/
  void
  OaFeasibilityChecker::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("Feasibility checker using OA cuts", RegisteredOptions::BonminCategory);
    roptions->AddStringOption2("feas_check_cut_types", "Choose the type of cuts generated when an integer feasible solution is found",
                               "outer-approx",
                               "outer-approx", "Generate a set of Outer Approximations cuts.",
                               "Benders", "Generate a single Benders cut.",
                               "If it seems too much memory is used should try Benders to use less");
    roptions->setOptionExtraInfo("feas_check_cut_types", 19);
    

    roptions->AddStringOption3("feas_check_discard_policy", "How cuts from feasibility checker are discarded",
                               "detect-cycles",
                               "detect-cycles", "Detect if a cycle occurs and only in this case force not to discard.",
                               "keep-all", "Force cuts from feasibility checker not to be discarded (memory hungry but sometimes better).",
                               "treated-as-normal", "Cuts from memory checker can be discarded as any other cuts (code may cycle then)",
                               "Normally to avoid cycle cuts from feasibility checker should not be discarded in the node where they are generated. "
                               "However Cbc sometimes does it if no care is taken which can lead to an infinite loop in Bonmin (usually on simple problems). "
                               "To avoid this one can instruct Cbc to never discard a cut but if we do that for all cuts it can lead to memory problems. "
                               "The default policy here is to detect cycles and only then impose to Cbc to keep the cut. "
                               "The two other alternative are to instruct Cbc to keep all cuts or to just ignore the problem and hope for the best");
    roptions->setOptionExtraInfo("feas_check_discard_policy", 19);

    roptions->AddLowerBoundedIntegerOption("generate_benders_after_so_many_oa", "Specify that after so many oa cuts have been generated Benders cuts should be generated instead.",
                                           0, 5000,
                                           "It seems that sometimes generating too many oa cuts slows down the optimization compared to Benders due to the size of the LP. "
                                           "With this option we specify that after so many OA cuts have been generated we should switch to Benders cuts.");
    roptions->setOptionExtraInfo("generate_benders_after_so_many_oa", 19);
  }

}/* End namespace Bonmin. */

// (C) Copyright International Business Machines (IBM) 2006, 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// P. Bonami, International Business Machines
//
// Date :  12/20/2006
//#define ECP_DEBUG
#include "BonEcpCuts.hpp"
#include "BonSolverHelp.hpp"
#include "BonBabInfos.hpp"
#include "BonCbc.hpp"
namespace Bonmin
{


  EcpCuts::EcpCuts(BabSetupBase & b):
      OaDecompositionBase(b,false, false)
  {
    assignLpInterface(NULL);
    b.options()->GetIntegerValue("ecp_max_rounds", numRounds_,b.prefix());
    b.options()->GetNumericValue("ecp_abs_tol", abs_violation_tol_,b.prefix());
    b.options()->GetNumericValue("ecp_rel_tol", rel_violation_tol_,b.prefix());
    b.options()->GetNumericValue("ecp_probability_factor", beta_,b.prefix());
  }

  double
  EcpCuts::doEcpRounds(OsiSolverInterface &si,
      bool leaveSiUnchanged,
      double* violation)
  {
    OsiSolverInterface * saveLp = lp_;
    lp_ = &si;
    OsiCuts cs;
    bool saveLeaveSi = leaveSiUnchanged_;
    leaveSiUnchanged_ = leaveSiUnchanged;
    generateCuts(si, cs);
    lp_ = saveLp;
    leaveSiUnchanged_ = saveLeaveSi;
    if (violation) *violation = violation_;
    return objValue_;
  }

  void
  EcpCuts::generateCuts(const OsiSolverInterface &si,
      OsiCuts & cs,
      const CglTreeInfo info) const
  {
    if (beta_ >=0) {
      BabInfo * babInfo = dynamic_cast<BabInfo *> (si.getAuxiliaryInfo());
      assert(babInfo);
      assert(babInfo->babPtr());
      const CbcNode * node = babInfo->babPtr()->model().currentNode();
      int level = (node == NULL) ? 0 : babInfo->babPtr()->model().currentNode()->depth();
      double rand = CoinDrand48();
      double score = pow(2.,-level)*beta_;
      //printf("depth %i, score %g , rand %g\n", level, score, rand);
      if (score <= rand)
        return;
    }
    double orig_violation = nlp_->getNonLinearitiesViolation(
          si.getColSolution(), si.getObjValue());
//#define ECP_DEBUG
#ifdef ECP_DEBUG
    std::cout<<"Initial Constraint violation: "<<orig_violation<<std::endl;
    std::cout<<"Initial objectvie value"<<si.getObjValue()<<std::endl;
#endif
    if (orig_violation <= abs_violation_tol_)
      return;
#ifdef ECP_DEBUG
    std::cout<<"Generating ECP cuts"<<std::endl;
#endif
    solverManip * lpManip = NULL;
    bool infeasible = false;
    violation_ = orig_violation;
    for (int i = 0 ; i < numRounds_ ; i++) {
      if ( violation_ > abs_violation_tol_ &&
          violation_ > rel_violation_tol_*orig_violation) {
        int numberCuts =  - cs.sizeRowCuts();
        const double * toCut = parameter().addOnlyViolated_ ?
            si.getColSolution():NULL;
        const OsiSolverInterface &localSi = (lpManip == NULL) ?
            si : *(lpManip->si());
        nlp_->getOuterApproximation(cs, localSi.getColSolution(), 1, toCut, parameter().global_);
        numberCuts += cs.sizeRowCuts();
        if (numberCuts > 0 && i + 1 < numRounds_ ) {
          if (lpManip==NULL) {
            assert(lp_ == NULL);
            if (lp_ == NULL)
              lpManip = new solverManip(si);
            else
              lpManip = new solverManip(lp_, true,true,
                  false,false);
          }
          installCuts(*lpManip->si(), cs,numberCuts);
#ifdef ECP_DEBUG
          std::cerr<<"Installed "<<numberCuts<<"cuts in lp"<<std::endl;
#endif
          lpManip->si()->resolve();
#ifdef ECP_DEBUG
          std::cerr<<"New objective "<<lpManip->si()->getObjValue()<<std::endl;
#endif
          if (lpManip->si()->isProvenPrimalInfeasible()) {
            infeasible = true;
#ifdef ECP_DEBUG
            std::cout<<"Stopping Ecp generation because problem became infeasible"<<std::endl;
#endif
            break;
          }
          violation_ =  nlp_->getNonLinearitiesViolation(
                lpManip->si()->getColSolution(),
                lpManip->si()->getObjValue());
#ifdef ECP_DEBUG
          std::cout<<"Constraint violation: "<<violation_<<std::endl;
#endif
        }
        else break;
      }
      else break;
    }
    if (!infeasible) {
      if (lpManip != NULL) {
        lpManip->si()->resolve();
        if (lpManip->si()->isProvenPrimalInfeasible())
          objValue_ = COIN_DBL_MAX;
        else
          objValue_ = lpManip->si()->getObjValue();
      }
    }
    else objValue_ = COIN_DBL_MAX;
    if (lpManip) {
      if (lp_ != NULL && lpManip != NULL) {
        lpManip->restore();
      }

      delete lpManip;
    }
#ifdef ECP_DEBUG
    std::cout<<"End ecp cut generation"<<std::endl;
#endif
    return;
  }

  void
  EcpCuts::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("ECP cuts generation", RegisteredOptions::BonminCategory);
    roptions->AddLowerBoundedIntegerOption("filmint_ecp_cuts",
        "Specify the frequency (in terms of nodes) at which some a la filmint ecp cuts are generated.",
        0,0,
        "A frequency of 0 amounts to to never solve the NLP relaxation.");
    roptions->setOptionExtraInfo("filmint_ecp_cuts",3);
    roptions->AddLowerBoundedIntegerOption
    ("ecp_max_rounds",
     "Set the maximal number of rounds of ECP cuts.",
     0,5,
     "");
    roptions->setOptionExtraInfo("ecp_max_rounds",3);
    roptions->AddLowerBoundedNumberOption
    ("ecp_abs_tol",
     "Set the absolute termination tolerance for ECP rounds.",
     0,false,1e-6,
     "");
    roptions->setOptionExtraInfo("ecp_abs_tol",3);
    roptions->AddLowerBoundedNumberOption
    ("ecp_rel_tol",
     "Set the relative termination tolerance for ECP rounds.",
     0,false,0.,
     "");
    roptions->setOptionExtraInfo("ecp_rel_tol",3);
    roptions->AddNumberOption
    ("ecp_probability_factor",
     "Factor appearing in formula for skipping ECP cuts.",
     10.,
     "Choosing -1 disables the skipping.");
    roptions->setOptionExtraInfo("ecp_probability_factor",3);
  }
} // end namespace bonmin.

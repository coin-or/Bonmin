
#include "BonLinearCutsGenerator.hpp"
#include "BonTMINLP2Quad.hpp"
#include "OsiClpSolverInterface.hpp"

//MILP cuts
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglTwomir.hpp"
#include "CglPreProcess.hpp"
#include "CglLandP.hpp"
#include "CglRedSplit.hpp"

namespace Bonmin {

  void LinearCutsGenerator::initialize(BabSetupBase & s){
    assert(dynamic_cast<TMINLP2TNLPQuadCuts *> (s.nonlinearSolver()->problem()));
    int freq;
    s.options()->GetIntegerValue("Gomory_cuts", freq,"bonmin.");
    if (freq) {
      Coin::SmartPtr<CuttingMethod> cg = new CuttingMethod;
      cg->frequency = freq;
      CglGomory * gom = new CglGomory;
      cg->cgl = gom;
      gom->setLimitAtRoot(512);
      gom->setLimit(50);
      cg->id = "Mixed Integer Gomory";
      methods_.push_back(cg);
    }

    s.options()->GetIntegerValue("mir_cuts",freq,"bonmin.");
    if (freq) {
      Coin::SmartPtr<CuttingMethod> cg = new CuttingMethod;
      cg->frequency = freq;
      CglMixedIntegerRounding2 * mir = new CglMixedIntegerRounding2;
      cg->cgl = mir;
      cg->id = "Mixed Integer Rounding";
      methods_.push_back(cg);


    }
    s.options()->GetIntegerValue("2mir_cuts",freq,"bonmin.");
    if (freq) {
      Coin::SmartPtr<CuttingMethod> cg = new CuttingMethod;
      cg->frequency = freq;
      CglTwomir * mir2 = new CglTwomir;
      cg->cgl = mir2;
      cg->id = "2-MIR";
      methods_.push_back(cg);
    }
    s.options()->GetIntegerValue("cover_cuts",freq,"bonmin.");
    if (freq) {
      Coin::SmartPtr<CuttingMethod> cg = new CuttingMethod;
      cg->frequency = freq;
      CglKnapsackCover * cover = new CglKnapsackCover;
      cg->cgl = cover;
      cg->id = "Cover";
      methods_.push_back(cg);
    }

    s.options()->GetIntegerValue("clique_cuts",freq,"bonmin.");
    if (freq) {
      Coin::SmartPtr<CuttingMethod> cg = new CuttingMethod;
      cg->frequency = freq;
      CglClique * clique = new CglClique;
      clique->setStarCliqueReport(false);
      clique->setRowCliqueReport(false);
      clique->setMinViolation(0.1);

      cg->cgl = clique;
      cg->id = "Clique";
      methods_.push_back(cg);
    }
    s.options()->GetIntegerValue("flow_cover_cuts",freq,"bonmin.");
    if (freq) {
      Coin::SmartPtr<CuttingMethod> cg = new CuttingMethod;
      cg->frequency = freq;
      CglFlowCover * flow = new CglFlowCover;
      cg->cgl = flow;
      cg->id = "Flow Covers";
      methods_.push_back(cg);
    }
    s.options()->GetIntegerValue("lift_and_project_cuts",freq,"bonmin.");
    if (freq) {
      Coin::SmartPtr<CuttingMethod> cg = new CuttingMethod;
      cg->frequency = freq;
      CglLandP * landp = new CglLandP;
      cg->cgl = landp;
      cg->id = "Lift-and-Project";
      methods_.push_back(cg);
    }
    s.options()->GetIntegerValue("reduce_and_split_cuts",freq,"bonmin.");
    if (freq) {
      Coin::SmartPtr<CuttingMethod> cg = new CuttingMethod;
      cg->frequency = freq;
      CglRedSplit * rands = new CglRedSplit;
      cg->cgl = rands;
      cg->id = "Reduce-and-Split";
      methods_.push_back(cg);
    }
  }

  void 
  LinearCutsGenerator::generateCuts(const OsiSolverInterface &solver, OsiCuts &cs,
                     const CglTreeInfo info) {

    //const OsiTMINLPInterface * tmp = dynamic_cast<const OsiTMINLPInterface *>(&solver);
    OsiTMINLPInterface * nlp = dynamic_cast<OsiTMINLPInterface *>(solver.clone());//const_cast<OsiTMINLPInterface *>(tmp);
    assert(nlp);
    OuterApprox oa;
    //si.writeMps("toto");
    int numberRows = nlp->getNumRows();
    for(int i = 0 ; i < 5 ; i++){
      nlp->resolve();
      OsiClpSolverInterface si;
      oa(*nlp, &si, solver.getColSolution(), true); 
      si.resolve();
      OsiCuts cuts;
      for(std::list<Coin::SmartPtr<CuttingMethod> >::const_iterator i = methods_.begin() ;
          i != methods_.end() ; i++){
         (*i)->cgl->generateCuts(si, cuts, info);
      }
      std::vector<OsiRowCut *> mycuts(cuts.sizeRowCuts());
      for(int i = 0 ; i < cuts.sizeRowCuts() ; i++){
        mycuts[i] = cuts.rowCutPtr(i);
        cs.insert(*mycuts[i]);
      }
      nlp->applyRowCuts((int)mycuts.size(), const_cast<const OsiRowCut **> (&mycuts[0]));
    }

    // Take off slack cuts
    std::vector<int> kept;
    int numberRowsNow = nlp->getNumRows();
    int * del = new int [numberRowsNow-numberRows];
    nlp->resolve();
    
    const double * activity = nlp->getRowActivity();
    const double * lb = nlp->getRowLower();
    const double * ub = nlp->getRowUpper();
    CoinRelFltEq eq(1e-06);
    //int nDelete=0;
    for (int i=numberRowsNow -1;i>=numberRows;i--) {
      if ( !(eq(activity[i], lb[i]) || eq(activity[i], ub[i])) )
        cs.eraseRowCut(i - numberRows);
    }
    delete [] del;
    delete nlp;
  }
}/* Ends Bonmin namespace.*/

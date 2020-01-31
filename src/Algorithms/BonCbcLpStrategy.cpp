// (C) Copyright Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 03/15/2006

#include "BonCbcLpStrategy.hpp"


// Cut generators
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglTwomir.hpp"
#include "CglPreProcess.hpp"
#include "CbcCutGenerator.hpp"

// Node selection
#include "CbcCompareActual.hpp"

#include "CbcBranchActual.hpp"


namespace Bonmin
{

  CbcStrategyChooseCuts::CbcStrategyChooseCuts():
    CbcStrategyDefault(),
    genFlag_(63)
  {
    CoinFillN(gen_freqs_,6,-99);
  }

  CbcStrategyChooseCuts::CbcStrategyChooseCuts(const CbcStrategyChooseCuts &other):
    CbcStrategyDefault(other),
    genFlag_(other.genFlag_)
  {
    CoinCopyN(other.gen_freqs_,6,gen_freqs_);
  }

  CbcStrategyChooseCuts::CbcStrategyChooseCuts(BabSetupBase &s,
                                               const std::string &prefix):
    CbcStrategyDefault(),
    genFlag_(0)
  {
    setup(s, prefix);
  }

  void CbcStrategyChooseCuts::setup(BabSetupBase &s,
                               const std::string &prefix){ 
    s.options()->GetIntegerValue("number_strong_branch", numberStrong_, prefix);
    s.options()->GetIntegerValue("number_before_trust", numberBeforeTrust_, prefix);

    int k = 0;
    
    bool set = s.options()->GetIntegerValue("probing_cuts", gen_freqs_[k], prefix);
    if(set==0) gen_freqs_[k] = -99;
    k++;

    set = s.options()->GetIntegerValue("Gomory_cuts", gen_freqs_[k], prefix);
    if(set==0) gen_freqs_[k] = -99;
    k++;
    
    set = s.options()->GetIntegerValue("cover_cuts", gen_freqs_[k], prefix);
    if(set==0) gen_freqs_[k] = -99;
    k++;
    
    set = s.options()->GetIntegerValue("clique_cuts", gen_freqs_[k], prefix);
    if(set==0) gen_freqs_[k] = -99;
    k++;
    
    set = s.options()->GetIntegerValue("flow_cover_cuts", gen_freqs_[k], prefix);
    if(set==0) gen_freqs_[k] = -99;
    k++;
    
    set = s.options()->GetIntegerValue("mir_cuts", gen_freqs_[k], prefix);
    if(set==0) gen_freqs_[k] = -99;
    k++;
    
  }

template<class X>
bool has_cg(CbcModel &model, const X& gen){
  int numberGenerators = model.numberCutGenerators();
  for (int iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
    CglCutGenerator * generator = model.cutGenerator(iGenerator)->generator();
    X * cgl = dynamic_cast<X *>(generator);
    if (cgl) {
      return true;
    }
  }
  return false;
}

#define ADD_CG(model, gen, setting, name) model.addCutGenerator(&gen,setting, name)

  void 
  CbcStrategyChooseCuts::setupCutGenerators(CbcModel &model){
    CglProbing probing;
    probing.setUsingObjective(true);
    probing.setMaxPass(1);
    probing.setMaxPassRoot(1);
    // Number of unsatisfied variables to look at
    probing.setMaxProbe(10);
    // How far to follow the consequences
    probing.setMaxLook(10);
    // Only look at rows with fewer than this number of elements
    probing.setMaxElements(200);
    probing.setMaxElementsRoot(300);
    //generator1.setRowCuts(3);

    CglGomory miG;
    // try larger limit
    miG.setLimit(300);

    CglKnapsackCover cover;

    CglClique clique;
    clique.setStarCliqueReport(false);
    clique.setRowCliqueReport(false);

    CglMixedIntegerRounding2 mixedGen;
    CglFlowCover flowGen;
    int k = 0;

    if(gen_freqs_[k]!= 0 && !has_cg(model, probing)){
      ADD_CG(model, probing, gen_freqs_[k], "Probing"); 
    }
    k++;

    if(gen_freqs_[k]!= 0 && !has_cg(model, miG)){
      ADD_CG(model, miG, gen_freqs_[k], "Gomory"); 
    }
    k++;
    
    if(gen_freqs_[k] != 0 && !has_cg(model, cover)){
      ADD_CG(model, cover, gen_freqs_[k], "Knapsack"); 
    }
    k++;

    if(gen_freqs_[k] != 0 && !has_cg(model, clique)){
      ADD_CG(model, clique, gen_freqs_[k], "Clique"); 
    }
    k++;

    if(gen_freqs_[k] != 0 && !has_cg(model, flowGen)){
      ADD_CG(model, flowGen, gen_freqs_[k], "FlowCover"); 
    }
    k++;

    if(gen_freqs_[k] != 0 && !has_cg(model, mixedGen)){
      ADD_CG(model, mixedGen, gen_freqs_[k], "MixedIntegerRounding2"); 
    }
  }

}

// (C) Copyright Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Common Public License.
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
#include "CglMixedIntegerRounding.hpp"
#include "CglTwomir.hpp"
#include "CglPreProcess.hpp"


// Node selection
#include "CbcCompareActual.hpp"

#include "CbcBranchActual.hpp"


namespace Bonmin
{
  CbcOaStrategy::CbcOaStrategy(int migFreq,
      int probFreq,
      int mirFreq,
      int coverFreq,
      int minReliability,
      int numberStrong,
      int nodeSelection,
      double intTol,
      int logLevel
                              ):
      CbcStrategy(),
      migFreq_(migFreq),
      probFreq_(probFreq),
      mirFreq_(mirFreq),
      coverFreq_(coverFreq),
      minReliability_(minReliability),
      numberStrong_(numberStrong),
      nodeSelection_(nodeSelection),
      intTol_(intTol),
      logLevel_(logLevel)
  {
    setPreProcessState(0);
  }

  CbcStrategy *
  CbcOaStrategy::clone () const
  {
    return new CbcOaStrategy(migFreq_, probFreq_,  mirFreq_, coverFreq_, minReliability_,
        numberStrong_, nodeSelection_, intTol_,
        logLevel_);
  }

  void
  CbcOaStrategy::setupCutGenerators(CbcModel & model)
  {

    CglGomory miGGen;

    CglProbing probGen;
    probGen.setUsingObjective(true);
    probGen.setMaxPass(3);
    probGen.setMaxProbe(100);
    probGen.setMaxLook(50);

    CglKnapsackCover knapsackGen;
    CglMixedIntegerRounding mixedGen;

    if (migFreq_ != 0)
      model.addCutGenerator(&miGGen,migFreq_,"GMI");
    if (probFreq_ != 0)
      model.addCutGenerator(&probGen,probFreq_,"Probing");
    if (coverFreq_ != 0)
      model.addCutGenerator(&knapsackGen,coverFreq_,"covers");
    if (mirFreq_ != 0)
      model.addCutGenerator(&mixedGen,mirFreq_,"MIR");

  }

/// Setup heuristics
  void
  CbcOaStrategy::setupHeuristics(CbcModel & model)
{}

/// Do printing stuff
  void
  CbcOaStrategy::setupPrinting(CbcModel & model,int modelLogLevel)
  {
    model.messageHandler()->setLogLevel(logLevel_);
    model.solver()->messageHandler()->setLogLevel(0);
    model.setPrintFrequency(100);
  }

// Other stuff e.g. strong branching
  void
  CbcOaStrategy::setupOther(CbcModel & model)
  {
    model.setNumberStrong(numberStrong_);
    model.setNumberBeforeTrust(minReliability_);

    model.setIntegerTolerance(intTol_);

    // Definition of node selection strategy
    CbcCompareObjective compare0;
    CbcCompareDepth compare1;
    CbcCompareDefault compare2;
    if (nodeSelection_==0) {
      model.setNodeComparison(compare0);
    }
    else if (nodeSelection_==1) {
      model.setNodeComparison(compare1);
    }
    else if (nodeSelection_==2) {
      compare2.setWeight(0.0);
      model.setNodeComparison(compare2);
    }
    else if (nodeSelection_==3) {
      model.setNodeComparison(compare2);
    }

  }

}

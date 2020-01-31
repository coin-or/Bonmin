// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date : 05/26/2005

#include "BonDummyHeuristic.hpp"
#include "CoinHelperFunctions.hpp"
#include "CbcModel.hpp"

#include "OsiAuxInfo.hpp"
namespace Bonmin
{
/// Default constructor
  DummyHeuristic::DummyHeuristic(CbcModel &model,
      OsiTMINLPInterface * si)
      :
      CbcHeuristic(model),
      nlp_(si)
  {}

  DummyHeuristic::DummyHeuristic(OsiTMINLPInterface * si)
      :
      CbcHeuristic(),
      nlp_(si)
  {}
/// Assign an OsiTMINLPInterface
  void
  DummyHeuristic::setNlp(OsiTMINLPInterface * si)
  {
    nlp_ = si;
  }
/// heuristic method
  int
  DummyHeuristic::solution(double &solutionValue, double *betterSolution)
  {
    OsiBabSolver * babSolver = dynamic_cast<OsiBabSolver *>
        (model_->solver()->getAuxiliaryInfo());
    //  double bestKnown = getObjValue();
    if (babSolver) {
      return babSolver->solution(solutionValue, betterSolution,
          model_->getNumCols());
    }
    return 0;
  }

}

// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date : 05/26/2005

#include "IpCbcDummyHeuristic.hpp"
#include "CoinHelperFunctions.hpp"
#include "CbcModel.hpp"

#include "OsiAuxInfo.hpp"
/// Default constructor
IpCbcDummyHeuristic::IpCbcDummyHeuristic(CbcModel &model,
    IpoptInterface * si)
    :
    CbcHeuristic(model),
    nlp_(si)
{}

IpCbcDummyHeuristic::IpCbcDummyHeuristic(IpoptInterface * si)
    :
    CbcHeuristic(),
    nlp_(si)
{}
/// Assign an IpoptInterface
void
IpCbcDummyHeuristic::assignInterface(IpoptInterface * si)
{
  nlp_ = si;
}
/// heuristic method
int
IpCbcDummyHeuristic::solution(double &solutionValue, double *betterSolution)
{
  OsiBabSolver * babSolver = dynamic_cast<OsiBabSolver *>
      (model_->solver()->getAuxiliaryInfo());
  //  double bestKnown = getObjValue();
  if(babSolver) {
    return babSolver->solution(solutionValue, betterSolution,
        model_->getNumCols());
  }
  return 0;
}


// (C) Copyright International Business Machines  2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Andreas Waechter          IBM       2007-09-01

#include "BonGuessHeuristic.hpp"
#include "CoinHelperFunctions.hpp"
#include "CbcModel.hpp"

//#include "OsiAuxInfo.hpp"
namespace Bonmin
{
  BonChooseVariable* GuessHeuristic::chooseMethod_ = NULL;
/// Default constructor
  GuessHeuristic::GuessHeuristic(CbcModel &model)
      :
    CbcHeuristic(model)
  {}

  /// heuristic method
  int
  GuessHeuristic::solution(double &solutionValue, double *betterSolution)
  {
    if (!chooseMethod_) {
      //printf("No choose method set yet\n");
      solutionValue = model_->getCurrentMinimizationObjValue();
      return -1;
    }
    int numberObjects = chooseMethod_->numberObjects();
    assert(numberObjects == model_->numberObjects());
    const double* upTotalChange = chooseMethod_->upTotalChange();
    const double* downTotalChange = chooseMethod_->downTotalChange();
    const int* upNumber = chooseMethod_->upNumber();
    const int* downNumber = chooseMethod_->downNumber();

    double sumUpTot = 0.;
    int numberUpTot = 0;
    double sumDownTot = 0.;
    int numberDownTot = 0;
    for (int i=0;i<numberObjects;i++) {
      sumUpTot += upTotalChange[i];
      numberUpTot += upNumber[i];
      sumDownTot += downTotalChange[i];
      numberDownTot += downNumber[i];
    }
    assert(numberUpTot && numberDownTot);
    double upAvrg=sumUpTot/numberUpTot;
    double downAvrg=sumDownTot/numberDownTot;

    OsiObject** object =  model_->objects();

    solutionValue = model_->getCurrentMinimizationObjValue();
    for (int iObj = 0; iObj < numberObjects; iObj++) {
      //printf("%3d upest=%e uptot=%e upnum=%d downest=%e downtot=%e downnum=%d  ", iObj, object[iObj]->upEstimate(), upTotalChange[iObj], upNumber[iObj], object[iObj]->downEstimate(), downTotalChange[iObj], downNumber[iObj]);
	
      double upEstimate = upNumber[iObj] ? object[iObj]->upEstimate()*upTotalChange[iObj]/upNumber[iObj] : object[iObj]->upEstimate()*upAvrg;
      double downEstimate = downNumber[iObj] ? object[iObj]->downEstimate()*downTotalChange[iObj]/downNumber[iObj] : object[iObj]->downEstimate()*downAvrg;
      //printf("up=%e down=%e\n", upEstimate, downEstimate);
      solutionValue += CoinMin(upEstimate,downEstimate);
    }
    //printf("solutionValue = %e\n", solutionValue);
    return -1;
  }

}

// Copyright (C) 2007, International Business Machines Corporation and others. 
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Joao P. Goncalves, International Business Machines Corporation
//
// Date : November 12, 2007

#include "CoinPragma.hpp"
#include "BonHeuristicDiveFractional.hpp"
#include "CbcModel.hpp"

namespace Bonmin
{
  HeuristicDiveFractional::HeuristicDiveFractional() 
    :
    HeuristicDive()
  {}

  HeuristicDiveFractional::HeuristicDiveFractional(BonminSetup * setup)
    :
    HeuristicDive(setup)
  {
    Initialize(setup->options());    
  }

  HeuristicDiveFractional::HeuristicDiveFractional(const HeuristicDiveFractional &copy)
    :
    HeuristicDive(copy)
  {}

  HeuristicDiveFractional & 
  HeuristicDiveFractional::operator=( const HeuristicDiveFractional& rhs)
  {
    if (this!=&rhs) {
      HeuristicDive::operator=(rhs);
    }
    return *this;
  }

  CbcHeuristic *
  HeuristicDiveFractional::clone() const
  {
    return new HeuristicDiveFractional(*this);
  }

  void
  HeuristicDiveFractional::setInternalVariables(TMINLP2TNLP* minlp)
  {
    // no variables to set
  }

  void
  HeuristicDiveFractional::selectVariableToBranch(TMINLP2TNLP* minlp,
						  const vector<int> & integerColumns,
						  const double* newSolution,
						  int& bestColumn,
						  int& bestRound)
  {
    double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
    
    const double* x_l = minlp->x_l();
    const double* x_u = minlp->x_u();

    // select a fractional variable to bound
    double smallestFraction = COIN_DBL_MAX;
    bestColumn = -1;
    bestRound = -1; // -1 rounds down, +1 rounds up
    for(int iIntCol=0; iIntCol<(int)integerColumns.size(); iIntCol++) {
      int iColumn = integerColumns[iIntCol];
      double value=newSolution[iColumn];
      if (fabs(floor(value+0.5)-value)>integerTolerance) {
	double below = floor(value);
	double downFraction = COIN_DBL_MAX;
	if(below >= x_l[iColumn])
	  downFraction = value-below;
	double above = ceil(value);
	double upFraction = COIN_DBL_MAX;
	if(above <= x_u[iColumn])
	  upFraction = ceil(value)-value;
	double fraction = 0;
	int round = 0;
	if(downFraction < upFraction) {
	  fraction = downFraction;
	  round = -1;
	} else if(downFraction > upFraction) {
	  fraction = upFraction;
	  round = 1;
	} else {
	  double randomNumber = CoinDrand48();
	  if(randomNumber<0.5) {
	    fraction = downFraction;
	    round = -1;
	  } else {
	    fraction = upFraction;
	    round = 1;
	  }	  
	}
	if(fraction<smallestFraction) {
	  smallestFraction = fraction;
	  bestColumn = iColumn;
	  bestRound = round;
	}
      }
    }

  }

  void
  HeuristicDiveFractional::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
    roptions->SetRegisteringCategory("Primal Heuristics", RegisteredOptions::BonminCategory);
   roptions->AddStringOption2(
     "heuristic_dive_fractional",
     "if yes runs the Dive Fractional heuristic",
     "no",
     "no", "",
     "yes", "",
     "");
   roptions->setOptionExtraInfo("heuristic_dive_fractional", 63);
  }

  void 
  HeuristicDiveFractional::Initialize(Ipopt::SmartPtr<Ipopt::OptionsList> options){
  }

}

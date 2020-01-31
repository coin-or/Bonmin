// Copyright (C) 2007, International Business Machines Corporation and others. 
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Joao P. Goncalves, International Business Machines Corporation
//
// Date : November 12, 2007

#include "CoinPragma.hpp"
#include "BonHeuristicDiveMIPFractional.hpp"
#include "CbcModel.hpp"

namespace Bonmin
{
#if 0
  HeuristicDiveMIPFractional::HeuristicDiveMIPFractional() 
    :
    HeuristicDiveMIP()
  {}
#endif

  HeuristicDiveMIPFractional::HeuristicDiveMIPFractional(BonminSetup * setup)
    :
    HeuristicDiveMIP(setup)
  {
    Initialize(setup->options());    
  }

  HeuristicDiveMIPFractional::HeuristicDiveMIPFractional(const HeuristicDiveMIPFractional &copy)
    :
    HeuristicDiveMIP(copy)
  {}

  HeuristicDiveMIPFractional & 
  HeuristicDiveMIPFractional::operator=( const HeuristicDiveMIPFractional& rhs)
  {
    if (this!=&rhs) {
      HeuristicDiveMIP::operator=(rhs);
    }
    return *this;
  }

  CbcHeuristic *
  HeuristicDiveMIPFractional::clone() const
  {
    return new HeuristicDiveMIPFractional(*this);
  }

  void
  HeuristicDiveMIPFractional::setInternalVariables(TMINLP2TNLP* minlp)
  {
    // no variables to set
  }

  void
  HeuristicDiveMIPFractional::selectVariableToBranch(TMINLP2TNLP* minlp,
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
  HeuristicDiveMIPFractional::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
    roptions->SetRegisteringCategory("Primal Heuristics", RegisteredOptions::BonminCategory);
   roptions->AddStringOption2(
     "heuristic_dive_MIP_fractional",
     "if yes runs the Dive MIP Fractional heuristic",
     "no",
     "no", "",
     "yes", "",
     "");
   roptions->setOptionExtraInfo("heuristic_dive_MIP_fractional", 63);
  }

  void 
  HeuristicDiveMIPFractional::Initialize(Ipopt::SmartPtr<Ipopt::OptionsList> options){
  }

}

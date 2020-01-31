// Copyright (C) 2007, International Business Machines Corporation and others. 
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Joao P. Goncalves, International Business Machines Corporation
//
// Date : November 12, 2007

#include "CoinPragma.hpp"
#include "BonHeuristicDiveVectorLength.hpp"
#include "CbcModel.hpp"

namespace Bonmin
{
  HeuristicDiveVectorLength::HeuristicDiveVectorLength() 
    :
    HeuristicDive(),
    columnLength_(NULL)
  {}

  HeuristicDiveVectorLength::HeuristicDiveVectorLength(BonminSetup * setup)
    :
    HeuristicDive(setup),
    columnLength_(NULL)
  {
    Initialize(setup->options());    
  }

  HeuristicDiveVectorLength::HeuristicDiveVectorLength(const HeuristicDiveVectorLength &copy)
    :
    HeuristicDive(copy),
    columnLength_(NULL)
  {
  }

  HeuristicDiveVectorLength & 
  HeuristicDiveVectorLength::operator=( const HeuristicDiveVectorLength& rhs)
  {
    if (this!=&rhs) {
      HeuristicDive::operator=(rhs);
      delete [] columnLength_;
      columnLength_ = NULL;
    }
    return *this;
  }

  CbcHeuristic *
  HeuristicDiveVectorLength::clone() const
  {
    return new HeuristicDiveVectorLength(*this);
  }

  void
  HeuristicDiveVectorLength::setInternalVariables(TMINLP2TNLP* minlp)
  {
    delete [] columnLength_;
    
    int numberColumns;
    int numberRows;
    int nnz_jac_g;
    int nnz_h_lag;
    Ipopt::TNLP::IndexStyleEnum index_style;
    minlp->get_nlp_info(numberColumns, numberRows, nnz_jac_g,
			nnz_h_lag, index_style);
    
    const double* x_sol = minlp->x_sol();

    // Get the indicies of the jacobian
    // This is also a way of knowing which variables are
    // used in each row
    int* indexRow = new int[nnz_jac_g];
    int* indexCol = new int[nnz_jac_g];
    minlp->eval_jac_g(numberColumns, x_sol, false,
		      numberRows, nnz_jac_g,
		      indexRow, indexCol, 0);
    columnLength_ = new int[numberColumns];
    int indexCorrection = (index_style == Ipopt::TNLP::C_STYLE) ? 0 : 1;
    int iniCol = -1;
    for(int i=0; i<nnz_jac_g; i++) {
      int thisIndexCol = indexCol[i]-indexCorrection;
      if(indexCol[i] != iniCol) {
	iniCol = indexCol[i];
	columnLength_[thisIndexCol] = 1;
      }
      else {
	columnLength_[thisIndexCol]++;
      }
    }
    
  }

  void
  HeuristicDiveVectorLength::selectVariableToBranch(TMINLP2TNLP* minlp,
						    const vector<int> & integerColumns,
						    const double* newSolution,
						    int& bestColumn,
						    int& bestRound)
  {
    double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);

    const double* x_l = minlp->x_l();
    const double* x_u = minlp->x_u();

    int numberColumns;
    int numberRows;
    int nnz_jac_g;
    int nnz_h_lag;
    Ipopt::TNLP::IndexStyleEnum index_style;
    minlp->get_nlp_info(numberColumns, numberRows, nnz_jac_g,
			nnz_h_lag, index_style);

    double* gradient_f = new double[numberColumns];

    double bestScore = COIN_DBL_MAX;
    bestColumn = -1;
    bestRound = -1; // -1 rounds down, +1 rounds up
    minlp->eval_grad_f(numberColumns,newSolution,true,gradient_f);
    for(int iIntCol=0; iIntCol<(int)integerColumns.size(); iIntCol++) {
      int iColumn = integerColumns[iIntCol];
      double value=newSolution[iColumn];
      if (fabs(floor(value+0.5)-value)>integerTolerance) {
	double below = floor(value);
	double downFraction = COIN_DBL_MAX;
	//double downCost = COIN_DBL_MAX;
	double gradient = gradient_f[iColumn];
	if(below >= x_l[iColumn])
	  downFraction = value-below;
	double above = ceil(value);
	double upFraction = COIN_DBL_MAX;
	if(above <= x_u[iColumn])
	  upFraction = ceil(value)-value;
	double objdelta = 0;
	int round = 0;
	if(gradient>=0.0 && upFraction < COIN_DBL_MAX) {
	  objdelta = gradient*upFraction;
	  round = 1;
	} else if(gradient<0.0 && downFraction < COIN_DBL_MAX) {
	  objdelta = gradient*downFraction;
	  round = -1;
	} else if(upFraction < COIN_DBL_MAX) {
	  objdelta = gradient*upFraction;
	  round = 1;
	} else {
	  objdelta = gradient*downFraction;
	  round = -1;
	}
	double score = (objdelta + 1e-6)/((double)columnLength_[iColumn]+1.0);
	if(score<bestScore) {
	  bestScore = score;
	  bestColumn = iColumn;
	  bestRound = round;
	}
      }
    }

    delete [] gradient_f;

  }

  void
  HeuristicDiveVectorLength::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
    roptions->SetRegisteringCategory("Primal Heuristics", RegisteredOptions::BonminCategory);
   roptions->AddStringOption2(
     "heuristic_dive_vectorLength",
     "if yes runs the Dive VectorLength heuristic",
     "no",
     "no", "",
     "yes", "",
     "");
   roptions->setOptionExtraInfo("heuristic_dive_vectorLength", 63);
  }

  void 
  HeuristicDiveVectorLength::Initialize(Ipopt::SmartPtr<Ipopt::OptionsList> options){
  }

}

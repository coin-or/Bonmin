// (C) Copyright International Business Machines Corporation 2007 
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Andreas Waechter, International Business Machines Corporation
//
// Date : 12/07/2007

#include "BonInitHeuristic.hpp"
#include "CoinHelperFunctions.hpp"

namespace Bonmin{
  
  InitHeuristic::InitHeuristic(double objValue, const double* sol,
			       CouenneProblem& cp):
    CbcHeuristic(),
    objValue_(COIN_DBL_MAX),
    sol_(NULL)
  {
    when_ = 1; // to be run at root

    setHeuristicName("InitHeuristic");
    nVars_ = cp.nVars();

    if (cp.checkNLP (sol, objValue)) {

      objValue_ = objValue;

      sol_ = new double[nVars_];

      CoinCopyN(sol, cp.nOrigVars(), sol_);
      cp.getAuxs(sol_);
    }
  }

  InitHeuristic::InitHeuristic(const InitHeuristic & other)
    :
    CbcHeuristic(other),
    objValue_(other.objValue_),
    nVars_(other.nVars_)
  {
    if (other.sol_) {
      sol_ = new double[nVars_];
      CoinCopyN(other.sol_, nVars_, sol_);
    }
    else {
      sol_ = NULL;
    }
  }
  
  CbcHeuristic * 
  InitHeuristic::clone() const{
    return new InitHeuristic(*this);
  }
  
  InitHeuristic &
  InitHeuristic::operator=(const InitHeuristic & rhs){
    if(this != &rhs){
      CbcHeuristic::operator=(rhs);
      objValue_ = rhs.objValue_;
      nVars_ = rhs.nVars_;
      if (sol_) {
	delete [] sol_;
	sol_ = NULL;
      }

      if (rhs.sol_) {
	sol_ = new double[nVars_];
	CoinCopyN(rhs.sol_, nVars_, sol_);
      }
    }
    return *this;
  }
  
  InitHeuristic::~InitHeuristic(){
    if(sol_)
      delete [] sol_;
  }
  
  int
  InitHeuristic::solution(double & objectiveValue, double * newSolution){

    if (!sol_) return 0;
    int retval = 0;
    if (objValue_ < objectiveValue) {
      CoinCopyN(sol_, nVars_, newSolution);
      objectiveValue = objValue_;
      retval = 1;
    }
    delete [] sol_;
    sol_ = NULL;
    return retval;
  }
}/** Ends namespace Bonmin.*/

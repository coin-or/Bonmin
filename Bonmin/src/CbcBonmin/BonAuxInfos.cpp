// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/23/2007

#include "BonAuxInfos.hpp"

namespace Bonmin{
  /** Default constructor.*/
  BabInfo::BabInfo(int type):
  OsiBabSolver(type),
  babPtr_(NULL),
  infeasibleNode_(false),
  nlpSolution_(NULL),
  numcols_(0),
  hasNlpSolution_(false){
  }
  
  /** Constructor from OsiBabSolver.*/
  BabInfo::BabInfo(const OsiBabSolver &other):
  OsiBabSolver(other),
  babPtr_(NULL),
  infeasibleNode_(false),
  nlpSolution_(NULL),
  numcols_(0),
  hasNlpSolution_(false){
  }
  
  /** Copy constructor.*/
  BabInfo::BabInfo(const BabInfo &other):
  OsiBabSolver(other),
  babPtr_(other.babPtr_),
  infeasibleNode_(other.infeasibleNode_),
  nlpSolution_(NULL),
  numcols_(other.numcols_),
  hasNlpSolution_(other.hasNlpSolution_){
    if(other.nlpSolution_!=NULL){
      assert(numcols_ > 0);
      nlpSolution_ = new double[numcols_ + 1];
      CoinCopyN(other.nlpSolution_, numcols_+1, nlpSolution_);
    } 
  }
  
  /** Destructor.*/
  BabInfo::~BabInfo(){
    if(nlpSolution_ != NULL)
      delete [] nlpSolution_;
  }
  
/** Virtual copy constructor.*/
OsiAuxInfo * 
  BabInfo::clone() const{
  return new BabInfo(*this);}
  
  /** Pass a solution found by an nlp solver.*/
  void 
  BabInfo::setNlpSolution(const double * sol, int numcols, double objValue){
    if(numcols_ < numcols){
      delete [] nlpSolution_;
      nlpSolution_ = NULL;}
    if(nlpSolution_ == NULL){
      nlpSolution_ = new double[numcols + 1];
      numcols_ = numcols;
    }
    CoinCopyN(sol,  numcols, nlpSolution_);
    nlpSolution_[numcols] = objValue;
  }
  
}/* end namespace Bonmin*/


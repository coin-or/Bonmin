// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/23/2007

#include "BonminConfig.h"
#include "BonAuxInfos.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"

namespace Bonmin
{
  /** Default constructor.*/
  AuxInfo::AuxInfo(int type):
      OsiBabSolver(type),
      infeasibleNode_(false),
      objValue_ (COIN_DBL_MAX),
      nlpSolution_(NULL),
      numcols_(0),
      hasNlpSolution_(false),
      bestSolution2_(make_referenced(std::vector<double>())),
      bestObj2_(make_referenced(COIN_DBL_MAX))
  {}

  /** Constructor from OsiBabSolver.*/
  AuxInfo::AuxInfo(const OsiBabSolver &other):
      OsiBabSolver(other),
      infeasibleNode_(false),
      objValue_ (COIN_DBL_MAX),
      nlpSolution_(NULL),
      numcols_(0),
      hasNlpSolution_(false),
      bestSolution2_(make_referenced(std::vector<double>())),
      bestObj2_(make_referenced(COIN_DBL_MAX))
  {}

  /** Copy constructor.*/
  AuxInfo::AuxInfo(const AuxInfo &other):
      OsiBabSolver(other),
      infeasibleNode_(other.infeasibleNode_),
      objValue_ (other.objValue_),
      nlpSolution_(NULL),
      numcols_(other.numcols_),
      hasNlpSolution_(other.hasNlpSolution_),
      bestSolution2_(other.bestSolution2_),
      bestObj2_(other.bestObj2_)
  {
    if (other.nlpSolution_!=NULL) {
      assert(numcols_ > 0);
      nlpSolution_ = new double[numcols_ + 1];
      CoinCopyN(other.nlpSolution_, numcols_+1, nlpSolution_);
    }
  }

  /** Destructor.*/
  AuxInfo::~AuxInfo()
  {
    if (nlpSolution_ != NULL)
      delete [] nlpSolution_;
  }

  /** Virtual copy constructor.*/
  OsiAuxInfo *
  AuxInfo::clone() const
  {
    return new AuxInfo(*this);
  }

  double AuxInfo::nlpObjValue ()
  {return hasNlpSolution_ ? objValue_ : COIN_DBL_MAX;}

  /** Pass a solution found by an nlp solver.*/
  void
  AuxInfo::setNlpSolution(const double * sol, int numcols, double objValue)
  {
    if (numcols_ < numcols) {
      delete [] nlpSolution_;
      nlpSolution_ = NULL;
    }
    if (nlpSolution_ == NULL) {
      nlpSolution_ = new double[numcols + 1];
      numcols_ = numcols;
    }
    CoinCopyN(sol,  numcols, nlpSolution_);
    nlpSolution_[numcols] = objValue;
    objValue_ = objValue;
  }

}/* end namespace Bonmin*/


// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/12/2007

#include "BonPseudoCosts.hpp"


namespace Bonmin
{

  PseudoCosts::PseudoCosts()
  {
    upTotalChange_ = NULL;
    downTotalChange_ = NULL;
    upNumber_ = NULL;
    downNumber_ = NULL;
    numberObjects_ = 0;
    numberBeforeTrusted_ = 0;
  }

  /** Copy constructor.*/
  PseudoCosts::PseudoCosts(const PseudoCosts & rhs)
  {
    upTotalChange_ = NULL;
    downTotalChange_ = NULL;
    upNumber_ = NULL;
    downNumber_ = NULL;
    numberObjects_ = rhs.numberObjects_;
    numberBeforeTrusted_ = rhs.numberBeforeTrusted_;
    if (numberObjects_ > 0) {
      upTotalChange_ = CoinCopyOfArray(rhs.upTotalChange_,numberObjects_);
      downTotalChange_ = CoinCopyOfArray(rhs.downTotalChange_,numberObjects_);
      upNumber_ = CoinCopyOfArray(rhs.upNumber_,numberObjects_);
      downNumber_ = CoinCopyOfArray(rhs.downNumber_,numberObjects_);
    }
  }


  /** Assignment operator const version.*/
  PseudoCosts &
  PseudoCosts::operator=(const PseudoCosts&rhs)
  {
    if (this != &rhs) {
      delete [] upTotalChange_;
      delete [] downTotalChange_;
      delete [] upNumber_;
      delete [] downNumber_;
      numberObjects_ = rhs.numberObjects_;
      numberBeforeTrusted_ = rhs.numberBeforeTrusted_;
      if (numberObjects_ > 0) {
        upTotalChange_ = CoinCopyOfArray(rhs.upTotalChange_,numberObjects_);
        downTotalChange_ = CoinCopyOfArray(rhs.downTotalChange_,numberObjects_);
        upNumber_ = CoinCopyOfArray(rhs.upNumber_,numberObjects_);
        downNumber_ = CoinCopyOfArray(rhs.downNumber_,numberObjects_);
      }
    }
    return *this;
  }

}/* End Bonmin namespace.*/


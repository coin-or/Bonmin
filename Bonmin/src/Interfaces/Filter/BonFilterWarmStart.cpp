// (C) Copyright International Business Machines Corporation, 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 11/21/2006

#include "BonFilterWarmStart.hpp"


namespace Bonmin{

  CoinWarmStartDiff *
  FilterWarmStart::generateDiff(const CoinWarmStart * const oldOne) const {

    const FilterWarmStart * old = dynamic_cast<const FilterWarmStart  *>(oldOne);

    if(size_ != old->size_){
      throw CoinError("Can not make difference for warm starts of differnet sizes",
		      "generateDiffs",
		      "FilterWarmStart");
    }
    FilterWarmStartDiff * diff = new FilterWarmStartDiff(size_);

    fint lastDiff=-1;
    for(fint i = 0 ; i < size_ ; i++){
      if(warmArray_[i] != old->warmArray_[i]){
	diff->differences.push_back(FilterWarmStartDiff::OneDiff(i, warmArray_[i] - old->warmArray_[i]));
	lastDiff = i;
      }
    }
    
    //    std::cout<<"Index of last difference "<<lastDiff<<", size "<<size_<<std::endl;
    diff->differences.resize(diff->differences.size());
    
    for(int i = 0 ; i < 14 ; i++){
      diff->istat_[i] = istat_[i];
    }
    return diff;
      
  }


  void
  FilterWarmStart::applyDiff(const CoinWarmStartDiff * diff){
    
    const FilterWarmStartDiff * diffF = dynamic_cast<const FilterWarmStartDiff  *>(diff);
    assert(diffF != NULL);
    fint end = diffF->differences.size();
    for(fint i = 0 ; i < end ; i++)
      {
	warmArray_[diffF->differences[i].first] += diffF->differences[i].second;
      }

    for(int i = 0 ; i < 14 ; i++)
      istat_[i] = diffF->istat_[i];
  }


FilterSolver::fint 
FilterWarmStart::def_istat[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};


} /* End namespace Bonmin */

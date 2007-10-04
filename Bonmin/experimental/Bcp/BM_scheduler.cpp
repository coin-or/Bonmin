// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 10/03/2007

#include "BM_scheduler.hpp"
#include <algorithm>
#include <cmath>

BM_scheduler::BM_scheduler():
    totalNumberIds_(0),
    freeIds_(),
    numFreeIds_(0),
    numNodeIds_(0){
}

BM_scheduler::BM_scheduler(const std::list<int> &freeIds):
    totalNumberIds_(0),
    freeIds_(freeIds),
    numFreeIds_(freeIds.size()),
    numNodeIds_(0){
}

void
BM_scheduler::add_free_ids(const std::list<int> &freeIds){
   int size = freeIds.size();
   totalNumberIds_ += size;
   numFreeIds_ += size;
   freeIds_.insert(freeIds_.begin(),freeIds.begin(), freeIds.end() );
}

int
BM_scheduler::request_sb_ids(int numIds, std::list<int> &ids){
   numIds = std::min(numIds, max_id_allocation());
   
   int numIdsToKeep = numFreeIds_ - numIds;
   while(numFreeIds_ > numIdsToKeep){
     ids.push_front(freeIds_.front());
     freeIds_.pop_front();
     numFreeIds_--;
   }
   return numIds;
}

int 
BM_scheduler::max_id_allocation(){
   double rho = 3;
  return floor(std::min((double) numFreeIds_, rho * (double) (totalNumberIds_ - numNodeIds_)/( (double) numNodeIds_) ));
}

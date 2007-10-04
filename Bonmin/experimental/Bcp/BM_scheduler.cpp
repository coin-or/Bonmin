// (C) 2007 Copyright International Business Machines Corporation
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
// Andreas Waechter, International Business Machines Corporation
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
BM_scheduler::setParams(double OverEstimationStatic,
			double SwitchToRateThreshold,
			double TimeRootNodeSolve,
			double FactorTimeHorizon,
			double OverEstimationRate){
  rho_static_ = OverEstimationStatic;
  switch_thresh_ = SwitchToRateThreshold;
  numSecRateInterval_ = (int)ceil(TimeRootNodeSolve*FactorTimeHorizon);
  request_counts_.reserve(numSecRateInterval_+1);
  std::fill_n(back_inserter(request_counts_), numSecRateInterval_+1, -1);
  request_counts_tot_ = 0;
  release_counts_.reserve(numSecRateInterval_+1);
  std::fill_n(back_inserter(release_counts_), numSecRateInterval_+1, -1);
  release_counts_tot_ = 0;
  counts_ptr_ = 0;
  time_last_action_ = 0;
  static_ = true;
  have_rates_ = false;
  rho_static_ = OverEstimationRate;
}

void
BM_scheduler::add_free_ids(const std::list<int> &freeIds){
  int size = freeIds.size();
  totalNumberIds_ += size;
  numFreeIds_ += size;
  freeIds_.insert(freeIds_.begin(),freeIds.begin(), freeIds.end() );
}

void
BM_scheduler::update_rates(int add_req, int add_rel)
{
  // Update the counts for the requests
  time_t time_now = time(NULL);
  if (time_now == time_last_action_) {
    request_counts_[counts_ptr_] += add_req;;
    release_counts_[counts_ptr_] += add_rel;;
  }
  else if (time_last_action_ == 0) {
    counts_ptr_ = 0;
    request_counts_[counts_ptr_] = add_req;
    release_counts_[counts_ptr_] = add_rel;
    time_last_action_ = time_now;
  }
  else {
    while (time_last_action_ < time_now) {
      request_counts_tot_ += request_counts_[counts_ptr_];
      release_counts_tot_ += release_counts_[counts_ptr_];
      counts_ptr_++;
      if (counts_ptr_ > numSecRateInterval_) {
	counts_ptr_ = 0;
	have_rates_ = true;
      }
      if (have_rates_) {
	request_counts_tot_ -= request_counts_[counts_ptr_];
	release_counts_tot_ -= release_counts_[counts_ptr_];
      }
      request_counts_[counts_ptr_] = 0;
      release_counts_[counts_ptr_] = 0;

      time_last_action_++;
    }
    request_counts_[counts_ptr_] = add_req;
    release_counts_[counts_ptr_] = add_rel;

    static_ = (!have_rates_ || numFreeIds_ >= (1.-switch_thresh_)*(totalNumberIds_-numNodeIds_));
  }
}

int
BM_scheduler::request_sb_ids(int numIds, std::list<int> &ids){

  // increase the count for requests by one
  update_rates(1, 0);
  
  numIds = std::min(numIds, max_id_allocation());
   
  int numIdsToKeep = numFreeIds_ - numIds;
  while(numFreeIds_ > numIdsToKeep){
    ids.push_front(freeIds_.front());
    freeIds_.pop_front();
    numFreeIds_--;
  }
  return numIds;
}

void
BM_scheduler::release_sb_id(int id){
  // increase the count for releases by one
  update_rates(0, 1);
  
  freeIds_.push_front(id);
  numFreeIds_++;
}

int 
BM_scheduler::max_id_allocation(){
  int retval;

  if (static_) {
    retval = (int)floor(std::min((double) numFreeIds_, rho_static_ * (double) (totalNumberIds_ - numNodeIds_)/( (double) numNodeIds_) ));
  }
  else {
    if (request_counts_tot_ == 0) {
      retval = numFreeIds_;
    }
    else {
      retval = std::min(numFreeIds_,(int)floor(rho_rate_*(double)(release_counts_tot_)/(double)(request_counts_tot_)));
    }
  }

  return retval;
}

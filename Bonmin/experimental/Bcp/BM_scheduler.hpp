// (C) 2007 Copyright International Business Machines Corporation
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
// Andreas Waechter, International Business Machines Corporation
//
// Date : 10/03/2007

#ifndef _BM_SCHEDULER_H
#define _BM_SCHEDULER_H
#include <list>
#include <vector>
#include "CoinTime.hpp"
#include <sys/time.h>

class BM_scheduler {
  /** Type of the container used to store the list of ids.*/
  typedef std::list<int> Id_Storage;

  /** Default constructor.*/
  BM_scheduler();

  /** Constructor with a list of free ids.*/
  BM_scheduler(const Id_Storage &freeIds);

  /** Method for setting scheduler parameters.
   * \param OverEstimationStatic: Factor for providing more IDs in static strategy.
   * \param SwitchToRateThreshold: When more than SwitchToRateThreshold times the number of strong-branching CPUs are busy, which to rate-based strategy.
   * \param TimeRootNodeSolve: Time to solve root node NLP (in secs)
   * \param FactorTimeHorizon: This number times TimeRootNodeSolve is used to compute the rates
   * \param OverEstimationRate: Factor for providing more IDs in rate-based strategy.
*/
  void setParams(double OverEstimationStatic,
		 double SwitchToRateThreshold,
		 double TimeRootNodeSolve,
		 double FactorTimeHorizon,
		 double OverEstimationRate);

  /** Pass in a list of freeIds_ to add.*/
  void add_free_ids(const Id_Storage &freeIds);

  /** Request for a number of id's to do some strong branching.
    * \param numIds : number of ids requested
    * \param ids : filled vector with the number of ids served.
    * \return number of ids served.
  */
  int request_sb_ids(int numIds, std::list<int> & ids);
  /** Request an id for processing nodes.
      \return id number or -1 if none is available. */
  inline int request_node_id(){
    if (numFreeIds_ == 0) return -1;
    numFreeIds_ --;
    numNodeIds_ ++;
    int id = freeIds_.front();
    freeIds_.pop_front();
    return id;
   }

  /** Gives back to scheduler an id used for strong branching.*/
  void release_sb_id(int id);

  /** Give back an id to scheduler used for processing a node
      (same as above but in addition decrement the number of free lp ids_)*/
  inline void release_node_id(int id){
    release_sb_id(id);
    numNodeIds_--;}

 private:
  /** Compute max allowed allocation of CPUs.*/
  int max_id_allocation();
  /** Update the counts and the static_ flag */
  void update_rates(int add_req, int add_rel);
 private:
  /** Store the total number of CPUs.*/
  int totalNumberIds_;
  /** List of free CPUs ids.*/
  Id_Storage freeIds_;
  /** number of freeIds in the list.*/
  int numFreeIds_;
  /** number of lp ids served.*/
  int numNodeIds_;
  /** overestimation factor for static strategy */
  double rho_static_;
  /** percentage threshold to swtich to rate-based strategy */
  double switch_thresh_;
  /** Number of seconds in time horizon for rate computation. */
  int numSecRateInterval_;
  /** vector for counting id requests per time unit */
  std::vector<int> request_counts_;
  /** total number of requests in considered time interval */
  int request_counts_tot_;
  /** vector for counting released sb id requests per time unit */
  std::vector<int> release_counts_;
  /** total number of releases in considered time interval */
  int release_counts_tot_;
  /** Array counter */
  int counts_ptr_;
  /** Time stamp of last request or release */
  time_t time_last_action_;
  /** overestimation factor for rate-based strategy */
  double rho_rate_;
  /** flag indicating whether we are in the static or the rate-based
   *  phase */
  bool static_;
  /** flag indicating whether we have rate information (i.e., the time
      horizon has passed at least once) */
  bool have_rates_;
};

#endif


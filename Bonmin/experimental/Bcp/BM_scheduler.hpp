// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 10/03/2007

#ifndef _BM_SCHEDULER_H
#define _BM_SCHEDULER_H
#include <list>
#include <vector>
#include "CoinTime.hpp"

class BM_scheduler {
  /** Type of the container used to store the list of ids.*/
  typedef std::list<int> Id_Storage;

  /** Default constructor.*/
  BM_scheduler();

  /** Constructor with a list of free ids.*/
  BM_scheduler(const Id_Storage &freeIds);

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
  inline void release_sb_id(int id){
    freeIds_.push_front(id);
    numFreeIds_++;}

  /** Give back an to scheduler id used for processing a node
      (same as above but in addition decrement the number of free lp ids_)*/
   inline void release_node_id(int id){
     release_sb_id(id);
     numNodeIds_--;}

 private:
   /** Compute max allowed allocation of CPUs.*/
   int max_id_allocation();
 private:
  /** Store the total number of CPUs.*/
  int totalNumberIds_;
  /** List of free CPUs ids.*/
  Id_Storage freeIds_;
  /** number of freeIds in the list.*/
  int numFreeIds_;
  /** number of lp ids served.*/
  int numNodeIds_;
};

#endif


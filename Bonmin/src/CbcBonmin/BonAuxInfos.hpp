// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/23/2007

#ifndef BonAuxInfos_H
#define BonAuxInfos_H
#include "OsiAuxInfo.hpp"
#include "BonCbc2.hpp"

namespace Bonmin {
class BabInfo : public OsiBabSolver {
public:
  BabInfo(int type):
  OsiBabSolver(type),
  babPtr_(NULL),
  infeasibleNode_(false){}
  
  BabInfo(const OsiBabSolver &other):
  OsiBabSolver(other){
  }
  BabInfo(const BabInfo &other):
  OsiBabSolver(other),
  babPtr_(other.babPtr_),
  infeasibleNode_(other.infeasibleNode_){
  }
  
  OsiAuxInfo * clone() const{
    return new BabInfo(*this);}
  void setBabPtr(Bab2 * babPtr){
    babPtr_ = babPtr;}
  
  Bab2 * babPtr(){
    return babPtr_;}
  
  void setFeasibleNode(){
    infeasibleNode_ = false;}
  
  void setInfeasibleNode(){
    infeasibleNode_ = true;}
  
  bool infeasibleNode(){
    return infeasibleNode_;}
protected: 
  /** Pointer to branch-and-bound algorithm.*/
  Bab2 * babPtr_;
  /** Say if current node was found infeasible during cut generation*/
  bool infeasibleNode_;
  
  };
}/* End namespace.*/

#endif


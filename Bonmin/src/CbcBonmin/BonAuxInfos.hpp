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
  babPtr_(NULL){}
  
  BabInfo(const OsiBabSolver &other):
  OsiBabSolver(other){
  }
  BabInfo(const BabInfo &other):
  OsiBabSolver(other),
  babPtr_(other.babPtr_){
  }
  
  OsiAuxInfo * clone() const{
    return new BabInfo(*this);}
  void setBabPtr(Bab2 * babPtr){
    babPtr_ = babPtr;}
  
  Bab2 * babPtr(){
    return babPtr_;}
protected: 
  Bab2 * babPtr_;
  };
}/* End namespace.*/

#endif


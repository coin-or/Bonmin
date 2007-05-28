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
#include <stdlib.h>
#include "OsiAuxInfo.hpp"
#include "BonCbc2.hpp"

namespace Bonmin {
  /** Bonmin class for passing info between components of branch-and-cuts.*/
class BabInfo : public OsiBabSolver {
public:
  /** Default constructor.*/
  BabInfo(int type);

  /** Constructor from OsiBabSolver.*/
  BabInfo(const OsiBabSolver &other);

  /** Copy constructor.*/
  BabInfo(const BabInfo &other);
  
  /** Destructor.*/
  virtual ~BabInfo();
  
  /** Virtual copy constructor.*/
  virtual OsiAuxInfo * clone() const;
  
  /** Set pointer to the branch-and-bound algorithm (to access CbcModel).*/
  void setBabPtr(Bab2 * babPtr){
    babPtr_ = babPtr;}
  
  /** Pointer to the branch-and-bound algorithm (to access CbcModel).*/
  Bab2 * babPtr(){
    return babPtr_;}
  
  /** Declare the node to be feasible.*/
  void setFeasibleNode(){
    infeasibleNode_ = false;}
  
  /** Declare the node to be infeasible.*/
  void setInfeasibleNode(){
    infeasibleNode_ = true;}
  
  /** Say if current node is found feasible by cut generators.*/
  bool infeasibleNode(){
    return infeasibleNode_;}
  
  /** Get solution found by nlp solver (or NULL if none found).*/
  const double * nlpSolution(){
    if(hasNlpSolution_)
      return nlpSolution_;
    else
      return NULL;
  }
    
  /** Pass a solution found by an nlp solver.*/
  void setNlpSolution(const double * sol, int numcols, double objValue);
  
  /** Say if has an nlp solution*/
  void setHasNlpSolution(bool b){
    hasNlpSolution_ = b;}
protected: 
  /** Pointer to branch-and-bound algorithm.*/
  Bab2 * babPtr_;
  /** Say if current node was found infeasible during cut generation*/
  bool infeasibleNode_;
  /** nlp solution found by heuristic if any.*/
  double * nlpSolution_;
  /** numcols_ gives the size of nlpSolution_.*/
  int numcols_;
  /** say if has a solution.*/
  bool hasNlpSolution_;
  
  };
}/* End namespace.*/

#endif


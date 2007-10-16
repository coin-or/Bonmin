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
#include <cstdlib>
#include <vector>
#include "OsiAuxInfo.hpp"
#include "CoinSmartPtr.hpp"


namespace Bonmin {

  //structure to store an object of class X in a Coin::ReferencedObject
  template<class X>
  struct SimpleReferenced : public Coin::ReferencedObject {
   /** The object.*/
   X object;
   const X& operator()() const{
     return object;}

   X& operator()() {
     return object;}

  };

   template <class X>
   SimpleReferenced<X>* make_referenced(const X& other){
     SimpleReferenced<X> * ret_val = new SimpleReferenced<X>;
     ret_val->object = other;
     return ret_val;
    }

  /** Bonmin class for passing info between components of branch-and-cuts.*/
class AuxInfo : public OsiBabSolver {
public:
  /** Default constructor.*/
  AuxInfo(int type);

  /** Constructor from OsiBabSolver.*/
  AuxInfo(const OsiBabSolver &other);

  /** Copy constructor.*/
  AuxInfo(const AuxInfo &other);
  
  /** Destructor.*/
  virtual ~AuxInfo();
  
  /** Virtual copy constructor.*/
  virtual OsiAuxInfo * clone() const;
  
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
  /** get the best solution computed with alternative objective function.*/
  const std::vector<double>& bestSolution2() const
  {
    return (*bestSolution2_)();
  }
  /** return objective value of the best solution computed with alternative
      objective function.*/
  double bestObj2() const
  {
    return (*bestObj2_)();
  }
  /** Set an alternate objective value.*/
  void setBestObj2(double o)
  {
    (*bestObj2_)() = o;
  }
  void setBestSolution2(int n, double * d)
  {
    (*bestSolution2_)().clear();
    (*bestSolution2_)().insert((*bestSolution2_)().end(),d, d+n);
  }
protected: 
  /** Say if current node was found infeasible during cut generation*/
  bool infeasibleNode_;
  /** nlp solution found by heuristic if any.*/
  double * nlpSolution_;
  /** numcols_ gives the size of nlpSolution_.*/
  int numcols_;
  /** say if has a solution.*/
  bool hasNlpSolution_;
  /** Stores the solution with alternate objective.*/
  Coin::SmartPtr< SimpleReferenced<std::vector<double> > > bestSolution2_;
  /** Alternate solution objective value.*/
  Coin::SmartPtr< SimpleReferenced<double> > bestObj2_;
  };
}/* End namespace.*/

#endif


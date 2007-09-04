// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 09/01/2007

#include "BonDiver.hpp"
#include "CoinFinite.hpp"
#include "CbcModel.hpp"
#define DIVE_DEBUG
namespace Bonmin {

/************************************************************************/
/*                CbcDiver methods                                       */
/************************************************************************/

      /// Default constructor.
  CbcDiver::CbcDiver(): CbcTree(),
			      nextOnBranch_(NULL){
  }

    ///Copy constructor.
  CbcDiver::CbcDiver(const CbcDiver &rhs):CbcTree(rhs),
				  nextOnBranch_(rhs.nextOnBranch_){
  }

    /// Assignment operator.
  CbcDiver & 
  CbcDiver::operator=(const CbcDiver &rhs){
    if(this != &rhs){
      CbcTree::operator=(rhs);
      nextOnBranch_ = rhs.nextOnBranch_;}
    return *this;
  }

    /// Destructor.
  CbcDiver::~CbcDiver(){
  }

    ///copy constructor.
  CbcTree * 
  CbcDiver::clone() const{
    return new CbcDiver(*this);}
    

    ///Return top node (next node to process.*/
  CbcNode * 
  CbcDiver::top() const{
#ifdef DIVE_DEBUG
    std::cerr<<"CbcDiver::top"<<std::endl;
#endif
    if(nextOnBranch_ != NULL){
      return nextOnBranch_;}
    else return CbcTree::top();
  }

  /// Add node to the heap.
  void 
  CbcDiver::push(CbcNode * x){
#ifdef DIVE_DEBUG
    std::cerr<<"CbcDiver::push"<<std::endl;
#endif
    if(x->branchingObject()->numberBranchesLeft() == x->branchingObject()->numberBranches()){//Not Backtracking
      assert(nextOnBranch_==NULL);//Should not happen twice in a row
      nextOnBranch_ = x;
    }
    else
      CbcTree::push(x);
  }

  /// Remove the top node of the heap.
  void 
  CbcDiver::pop(){
#ifdef DIVE_DEBUG
    std::cerr<<"CbcDiver::pop"<<std::endl;
#endif
    if(nextOnBranch_ != NULL){
      nextOnBranch_ = NULL;}
    else
      CbcTree::pop();
  }
  /// Remove the best node from the heap and return it
  CbcNode * 
  CbcDiver::bestNode(double cutoff){
#ifdef DIVE_DEBUG
    std::cerr<<"CbcDiver::bestNode"<<std::endl;
#endif
    if(nextOnBranch_ != NULL){
      if(nextOnBranch_->objectiveValue() < cutoff){
      CbcNode * ret_val = nextOnBranch_;
      nextOnBranch_ = NULL;
      return ret_val;
      }
      else nextOnBranch_ = NULL;
    }
    return CbcTree::bestNode(cutoff);
  }


    /** Test if empty. */
  bool CbcDiver::empty(){
    return (CbcTree::empty() && (nextOnBranch_ == NULL) );
  }

    /*! \brief Prune the tree using an objective function cutoff
      if nextOnBranch_ exists we push it on the heap and call CbcTree function
    */
    void 
    CbcDiver::cleanTree(CbcModel * model, double cutoff, double & bestPossibleObjective)
    {
      if(nextOnBranch_ != NULL){
	CbcTree::push(nextOnBranch_);
	nextOnBranch_=NULL;
      }
      CbcTree::cleanTree(model,cutoff, bestPossibleObjective);
    }

    /// Get best possible objective function in the tree
    double 
    CbcDiver::getBestPossibleObjective(){
      double bestPossibleObjective = (nextOnBranch_ != NULL) ? nextOnBranch_->objectiveValue() : 1e100;
      for(unsigned int i = 0 ; i < nodes_.size() ; i++){
	if(nodes_[i] == NULL) continue;
	const double & obj = nodes_[i]->objectiveValue();
	if(obj < bestPossibleObjective){
	  bestPossibleObjective = obj;
	}
      }
      return bestPossibleObjective;
    }


/************************************************************************/
/*                CbcDfsDiver methods                                    */
/************************************************************************/
      /// Default constructor.
  CbcDfsDiver::CbcDfsDiver(): CbcTree(),
			      dive_(),
                              diveListSize_(0),
                              divingBoardDepth_(-1),
                              cutoff_(1e100),
                              nBacktracks_(0),
                              maxDiveBacktracks_(5),
                              maxDiveDepth_(COIN_INT_MAX){
  }

    ///Copy constructor.
  CbcDfsDiver::CbcDfsDiver(const CbcDfsDiver &rhs):CbcTree(rhs),
			                           dive_(rhs.dive_),
                                                   diveListSize_(rhs.diveListSize_),
                                                   divingBoardDepth_(rhs.divingBoardDepth_),
                                                   cutoff_(rhs.cutoff_),
                                                   nBacktracks_(rhs.nBacktracks_),
                                                   maxDiveBacktracks_(rhs.maxDiveBacktracks_),
                                                   maxDiveDepth_(maxDiveDepth_){
  }

    /// Assignment operator.
  CbcDfsDiver & 
  CbcDfsDiver::operator=(const CbcDfsDiver &rhs){
    if(this != &rhs){
      CbcTree::operator=(rhs);
      dive_ = rhs.dive_;
      diveListSize_ = rhs.diveListSize_;
      divingBoardDepth_ = rhs.divingBoardDepth_;
      cutoff_ = rhs.cutoff_;
      nBacktracks_ = rhs.nBacktracks_;
      maxDiveBacktracks_ = rhs.maxDiveBacktracks_;
      maxDiveDepth_ = maxDiveDepth_;}
    return *this;
  }

    /// Destructor.
  CbcDfsDiver::~CbcDfsDiver(){
  }

    ///copy constructor.
  CbcTree * 
  CbcDfsDiver::clone() const{
    return new CbcDfsDiver(*this);}

  ///Return top node (next node to process.*/
  CbcNode * 
  CbcDfsDiver::top() const{
#ifdef DIVE_DEBUG
    std::cerr<<"CbcDfsDiver::top"<<std::endl;
#endif
    if(mode_ != CbcDfsDiver::FindSolutions){
      assert(dive_.empty());
      CbcTree::top();
    }
    if(diveListSize_){
      return dive_.front();}
    else return CbcTree::top();
  }

  /// Add node to the heap.
  void 
  CbcDfsDiver::push(CbcNode * x){
#ifdef DIVE_DEBUG
    std::cerr<<"CbcDfsDiver::push"<<std::endl;
#endif
    if(mode_ != CbcDfsDiver::FindSolutions){
      CbcTree::push(x);
      assert(dive_.empty());
    }
    //Always push on dive;
    dive_.push_front(x);
    diveListSize_++;
  }

  /// Remove the top node of the heap.
  void 
  CbcDfsDiver::pop(){
#ifdef DIVE_DEBUG
    std::cerr<<"CbcDfsDiver::pop"<<std::endl;
#endif
    if(mode_ != CbcDfsDiver::FindSolutions){
      assert(dive_.empty());
    }
    if(!dive_.empty()){
      dive_.pop_front();
      diveListSize_--;}
    else
      CbcTree::pop();
  }

  /// Remove the best node from the heap and return it
  CbcNode * 
  CbcDfsDiver::bestNode(double cutoff){
    if(mode_ != CbcDfsDiver::FindSolutions){
      assert(dive_.empty());
      CbcTree::bestNode(cutoff);
    }
#ifdef DIVE_DEBUG
    std::cerr<<"CbcDfsDiver::bestNode"<<std::endl;
#endif
    assert(nBacktracks_ < maxDiveBacktracks_);
    CbcNode * node = NULL;
    while(diveListSize_ > 0){
       node = dive_.front();
       dive_.pop_front();
       diveListSize_ --;
       assert(node);
       assert((node->depth() - divingBoardDepth_) <= maxDiveDepth_);
       if(node->objectiveValue() > cutoff){//throw away node for now just put it on the heap as deleting a node is
                                           //more complicated than that (has to delete nodeInfo, cuts...)
#ifdef DIVE_DEBUG
         std::cerr<<"CbcDfsDiver::bestNode"
                  <<", node above cutoff"<<std::endl;
#endif
         CbcTree::push(node);
         node = NULL;
         nBacktracks_++;
       }
       else if(node->guessedObjectiveValue() > cutoff){//Put it on the real heap
#ifdef DIVE_DEBUG
         std::cerr<<"CbcDfsDiver::bestNode"
                  <<", node estimates above cutoff"<<std::endl;
#endif
         CbcTree::push(node);
         nBacktracks_++;
         node = NULL;
       }
       else if((node->depth() - divingBoardDepth_) > maxDiveDepth_){//Put it on the real heap
#ifdef DIVE_DEBUG
         std::cerr<<"CbcDfsDiver::bestNode"
                  <<", node too deep"<<std::endl;
#endif
         CbcTree::push(node);
         nBacktracks_++;
         node = NULL;
       }
      else if(node->branchingObject()->numberBranchesLeft() < node->branchingObject()->numberBranches()){//Backtracking
        nBacktracks_++;
#ifdef DIVE_DEBUG
         std::cerr<<"CbcDfsDiver::bestNode"
                  <<", backtracking"<<std::endl;
#endif
      }
      if(nBacktracks_ >= maxDiveBacktracks_){//Push all the node in dive_ onto nodes_
#ifdef DIVE_DEBUG
         std::cerr<<"CbcDfsDiver::bestNode"
                  <<", maximum number of backtracks attained emptying dive_"<<std::endl;
#endif
        pushDiveOntoHeap(-COIN_DBL_MAX);
        if(node != NULL) CbcTree::push(node);
        node = NULL;
      }
      if(node != NULL)
        return node;
    }
    assert(node == NULL);
    assert(dive_.empty());
    assert(diveListSize_ == 0);
    node = CbcTree::bestNode(cutoff);
    divingBoardDepth_ = node->depth();
    nBacktracks_ = 0;
    return node;
  }


  bool CbcDfsDiver::pushDiveOntoHeap(double cutoff){
    while(!dive_.empty() && (dive_.front() == NULL || dive_.front()->objectiveValue() > cutoff)){
        if(dive_.front() != NULL) CbcTree::push(dive_.front());
        dive_.pop_front();
        diveListSize_--;
    }
  }
  /** Test if empty. */
  bool CbcDfsDiver::empty(){
    return (CbcTree::empty() && dive_.empty());
  }

    /*! \brief Prune the tree using an objective function cutoff
      if nextOnBranch_ exists we push it on the heap and call CbcTree function
    */
    void 
    CbcDfsDiver::cleanTree(CbcModel * model, double cutoff, double & bestPossibleObjective)
    {
      pushDiveOntoHeap(cutoff);
      CbcTree::cleanTree(model,cutoff, bestPossibleObjective);
    }

    /// Get best possible objective function in the tree
    double 
    CbcDfsDiver::getBestPossibleObjective(){
#ifdef DIVE_DEBUG
         std::cerr<<"CbcDfsDiver::getBestPossibleObjective"<<std::endl;
#endif
      double bestPossibleObjective = CbcTree::empty() ? COIN_DBL_MAX : CbcTree::getBestPossibleObjective();
      for(std::list<CbcNode *>::iterator i = dive_.begin() ; i != dive_.end() ; i++){
	if(*i == NULL) continue;
	const double & obj = (*i)->objectiveValue();
	if(obj < bestPossibleObjective){
	  bestPossibleObjective = obj;
	}
      }
      return bestPossibleObjective;
    }
#ifdef COIN_HAS_BONMIN
    ///Register the options of the method.
    void 
    CbcDfsDiver::registerOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions){
      roptions->AddLowerBoundedIntegerOption("max-backtracks-in-dive",
                                                  "Set the number of backtracks in a dive when using dfs-dive tree search"
                                                  "strategy.",
                                                  0,5,
                                                  "");

      roptions->AddLowerBoundedIntegerOption("max-dive-depth",
                                                   "When using dfs-dive search. Maximum depth to go to from the diving "
                                                   "board (node where the diving started.",
                                                   0,INT_MAX,
                                                   "");

    }
   
    /// Initialize the method (get options)
    void 
    CbcDfsDiver::initialize(Ipopt::SmartPtr<Ipopt::OptionsList> options){
      options->GetIntegerValue("max-dive-depth", maxDiveDepth_,"bonmin.");
      options->GetIntegerValue("max-backtracks-in-dive", maxDiveBacktracks_,"bonmin.");
    }

  /** Changes the mode of comparison of the tree for "safety reasons" if the mode really changes we always 
      finish the current dive and put all the node back onto the heap.*/
  void 
  CbcDfsDiver::setComparisonMode(ComparisonModes newMode){
    if(newMode != mode_){
       mode_ = newMode;
       //Empty heap
       pushDiveOntoHeap(-COIN_DBL_MAX);
    }
  }



  // This allows any method to change behavior as it is called
  // after each solution
  void 
  DiverCompare::newSolution(CbcModel * model){
    assert(diver_);
    if(model->getSolutionCount() > numberSolToStopDive_ && diver_->getComparisonMode() != CbcDfsDiver::FindSolutions){
      diver_->setComparisonMode(CbcDfsDiver::CloseBound);
    }
  }

  /// This is test function
  bool 
  DiverCompare::test (CbcNode * x, CbcNode * y){
    assert(diver_);
    assert(comparisonDive_);
    assert(comparisonBound_);
    CbcDfsDiver::ComparisonModes mode = diver_->getComparisonMode();
    if(mode == CbcDfsDiver::FindSolutions){
      comparisonDive_->test(x,y);}
    else if(mode == CbcDfsDiver::CloseBound){
      comparisonBound_->test(x,y);}
    else if(mode == CbcDfsDiver::LimitTreeSize){
       comparisonDepth_.test(x,y);}
  }


  // This Also allows any method to change behavior as it is called
  // after each solution
  void 
  DiverCompare::newSolution(CbcModel * model,
			   double objectiveAtContinuous,
			   int numberInfeasibilitiesAtContinuous){
    newSolution(model);}

  // This allows any method to change behavior as it is called
  // after every 1000 nodes.
  // Return true if want tree re-sorted
  bool 
  DiverCompare::every1000Nodes(CbcModel * model,int numberNodes){
    assert(diver_);
    if(numberNodes > numberNodesToLimitTreeSize_  && diver_->getComparisonMode() != CbcDfsDiver::LimitTreeSize){
      diver_->setComparisonMode(CbcDfsDiver::LimitTreeSize);
      return true;
    }
  }

}/* Ends namespace Bonmin.*/
#endif


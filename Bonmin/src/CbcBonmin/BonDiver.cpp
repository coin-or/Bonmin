// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 09/01/2007

#include "BonDiver.hpp"
#include "CoinFinite.hpp"
#include "CbcModel.hpp"
#include "CbcConfig.h"
#include "BonBabSetupBase.hpp"
//#define DIVE_DEBUG
namespace Bonmin
{

  /************************************************************************/
  /*                CbcDiver methods                                       */
  /************************************************************************/

  /// Default constructor.
  CbcDiver::CbcDiver(): CbcTree(),
      treeCleaning_(false),
      nextOnBranch_(NULL),
      stop_diving_on_cutoff_(false)
  {}

  ///Copy constructor.
  CbcDiver::CbcDiver(const CbcDiver &rhs):CbcTree(rhs),
      treeCleaning_(rhs.treeCleaning_),
      nextOnBranch_(rhs.nextOnBranch_),
      stop_diving_on_cutoff_(rhs.stop_diving_on_cutoff_)
  {}

  /// Assignment operator.
  CbcDiver &
  CbcDiver::operator=(const CbcDiver &rhs)
  {
    if (this != &rhs) {
      CbcTree::operator=(rhs);
      treeCleaning_ = rhs.treeCleaning_;
      nextOnBranch_ = rhs.nextOnBranch_;
      stop_diving_on_cutoff_ = rhs.stop_diving_on_cutoff_;
    }
    return *this;
  }

  /// Destructor.
  CbcDiver::~CbcDiver()
  {}

  ///copy constructor.
  CbcTree *
  CbcDiver::clone() const
  {
    return new CbcDiver(*this);
  }


  ///Return top node (next node to process.*/
  CbcNode *
  CbcDiver::top() const
  {
#ifdef DIVE_DEBUG
    std::cout<<"CbcDiver::top"<<std::endl;
#endif
    if (nextOnBranch_ != NULL && !treeCleaning_) {
      return nextOnBranch_;
    }
    else return CbcTree::top();
  }

  /// Add node to the heap.
  void
  CbcDiver::push(CbcNode * x)
  {
#ifdef DIVE_DEBUG
    std::cout<<"CbcDiver::push"<<std::endl;
#endif
    if (treeCleaning_) return CbcTree::push(x);

    if (x->branchingObject()->numberBranchesLeft() == x->branchingObject()->numberBranches()) {//Not Backtracking
      assert(nextOnBranch_==NULL);//Should not happen twice in a row
      nextOnBranch_ = x;
    }
    else
      CbcTree::push(x);
  }

  /// Remove the top node of the heap.
  void
  CbcDiver::pop()
  {
#ifdef DIVE_DEBUG
    std::cout<<"CbcDiver::pop"<<std::endl;
#endif
    if (nextOnBranch_ != NULL && !treeCleaning_) {
      nextOnBranch_ = NULL;
    }
    else
      CbcTree::pop();
  }
  /// Remove the best node from the heap and return it
  CbcNode *
  CbcDiver::bestNode(double cutoff)
  {
#ifdef DIVE_DEBUG
    std::cout<<"CbcDiver::bestNode"<<std::endl;
#endif
    if (nextOnBranch_ != NULL && !treeCleaning_) {
      if (nextOnBranch_->objectiveValue() < cutoff) {
        if (stop_diving_on_cutoff_ && nextOnBranch_->guessedObjectiveValue() >= cutoff) {
          //printf("Stopping dive %g %g\n", nextOnBranch_->objectiveValue(), nextOnBranch_->guessedObjectiveValue());
          CbcTree::push(nextOnBranch_);
          nextOnBranch_ = NULL;
          return CbcTree::bestNode(cutoff);
        }
        //printf("Diving!! %g %g\n", nextOnBranch_->objectiveValue(), nextOnBranch_->guessedObjectiveValue());
        CbcNode * ret_val = nextOnBranch_;
        nextOnBranch_ = NULL;
        return ret_val;
      }
      assert(true && "We did not think we get here.");
      CbcTree::push(nextOnBranch_);//For safety
      nextOnBranch_ = NULL;
    }
    return CbcTree::bestNode(cutoff);
  }


  /** Test if empty. */
  bool CbcDiver::empty()
  {
    return (CbcTree::empty() && (nextOnBranch_ == NULL) );
  }

  /*! \brief Prune the tree using an objective function cutoff
    if nextOnBranch_ exists we push it on the heap and call CbcTree function
  */
  void
  CbcDiver::cleanTree(CbcModel * model, double cutoff, double & bestPossibleObjective)
  {
    if (nextOnBranch_ != NULL)
      CbcTree::push(nextOnBranch_);
    nextOnBranch_=NULL;
    treeCleaning_ = true;
    CbcTree::cleanTree(model,cutoff, bestPossibleObjective);
    treeCleaning_ = false;
  }

  /// Get best possible objective function in the tree
  double
  CbcDiver::getBestPossibleObjective()
  {
    double bestPossibleObjective = (nextOnBranch_ != NULL) ? nextOnBranch_->objectiveValue() : 1e100;
    for (unsigned int i = 0 ; i < nodes_.size() ; i++) {
      if (nodes_[i] == NULL) continue;
      const double & obj = nodes_[i]->objectiveValue();
      if (obj < bestPossibleObjective) {
        bestPossibleObjective = obj;
      }
    }
    return bestPossibleObjective;
  }

  void
  CbcDiver::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("Diving options", RegisteredOptions::UndocumentedCategory);
    roptions->AddStringOption2(
      "stop_diving_on_cutoff",
      "Flag indicating whether we stop diving based on guessed feasible objective and the current cutoff",
      "no",
      "no", "",
      "yes", "");
    roptions->setOptionExtraInfo("stop_diving_on_cutoff", 63);

  }

  /// Initialize the method (get options)
  void
  CbcDiver::initialize(BabSetupBase &b)
  {
    b.options()->GetBoolValue("stop_diving_on_cutoff", stop_diving_on_cutoff_,
        b.prefix());
  }


  /************************************************************************/
  /*                CbcProbedDiver methods                                       */
  /************************************************************************/

  /// Default constructor.
  CbcProbedDiver::CbcProbedDiver(): CbcTree(),
      treeCleaning_(false),
      nextOnBranch_(NULL),
      candidateChild_(NULL),
      stop_diving_on_cutoff_(false)
  {}

  ///Copy constructor.
  CbcProbedDiver::CbcProbedDiver(const CbcProbedDiver &rhs):CbcTree(rhs),
      treeCleaning_(rhs.treeCleaning_),
      nextOnBranch_(rhs.nextOnBranch_),
      candidateChild_(rhs.candidateChild_),
      stop_diving_on_cutoff_(rhs.stop_diving_on_cutoff_)
  {}

  /// Assignment operator.
  CbcProbedDiver &
  CbcProbedDiver::operator=(const CbcProbedDiver &rhs)
  {
    if (this != &rhs) {
      CbcTree::operator=(rhs);
      treeCleaning_ = rhs.treeCleaning_;
      nextOnBranch_ = rhs.nextOnBranch_;
      candidateChild_ = rhs.candidateChild_;
      stop_diving_on_cutoff_ = rhs.stop_diving_on_cutoff_;
    }
    return *this;
  }

  /// Destructor.
  CbcProbedDiver::~CbcProbedDiver()
  {}

  ///copy constructor.
  CbcTree *
  CbcProbedDiver::clone() const
  {
    return new CbcProbedDiver(*this);
  }


  ///Return top node (next node to process.*/
  CbcNode *
  CbcProbedDiver::top() const
  {
#ifdef DIVE_DEBUG
    std::cout<<"CbcProbedDiver::top"<<std::endl;
#endif
    if (nextOnBranch_ != NULL && !treeCleaning_) {
      return nextOnBranch_;
    }
    else if(candidateChild_ != NULL && !treeCleaning_){
      return candidateChild_;
    }
    else return CbcTree::top();
  }

  /// Add node to the heap.
  void
  CbcProbedDiver::push(CbcNode * x)
  {
#ifdef DIVE_DEBUG
    std::cout<<"CbcProbedDiver::push"<<std::endl;
#endif
    if (treeCleaning_) return CbcTree::push(x);

    if (x->branchingObject()->numberBranchesLeft() == x->branchingObject()->numberBranches()) {
      //Not Backtracking
      if(nextOnBranch_ == NULL && candidateChild_ == NULL){
#ifdef DIVE_DEBUG
        printf("Putting parent %i. objective value %g\n", x, x->objectiveValue());
#endif
        nextOnBranch_ = x;
        return;
      }
      if(candidateChild_ == NULL){
#ifdef DIVE_DEBUG
        printf("Putting first kid %i. objective value %g\n", x, x->objectiveValue());
#endif
        candidateChild_ = x;
        return;}
#ifdef DIVE_DEBUG
        printf("Putting second kid %i. objective value %g\n", x, x->objectiveValue());
#endif
      // If we get here we have two nodes follow on the best one and put the other on the tree
      if(comparison_.compareNodes(x,candidateChild_)){// Follow on candidateChild_
#ifdef DIVE_DEBUG
        printf("first child %i is found better.\n", x);
#endif
        nextOnBranch_ = candidateChild_;
        CbcTree::push(x);
        candidateChild_ = NULL;
        return;
      }
      else {// Follow on x
#ifdef DIVE_DEBUG
        printf("second child %i is found better\n", x);
#endif
        nextOnBranch_ = x;
        CbcTree::push(candidateChild_);
        candidateChild_ = NULL;
        return;
      }
    }
    else {
#ifdef DIVE_DEBUG
        printf("Putting back parent %i.\n",x);
#endif
      if(nextOnBranch_ != NULL){
#ifdef DIVE_DEBUG
        printf("Putting nextOnBranch_ in candidateChild_.\n");
#endif
        assert(candidateChild_ == NULL);
        candidateChild_ = nextOnBranch_;
      }
      nextOnBranch_ = x;
      return;
      //CbcTree::push(x);
    }
  }
  /// Remove the top node of the heap.
  void
  CbcProbedDiver::pop()
  {
#ifdef DIVE_DEBUG
    std::cout<<"CbcProbedDiver::pop"<<std::endl;
#endif
    if (nextOnBranch_ != NULL && !treeCleaning_) {
      nextOnBranch_ = NULL;
    }
    else if(candidateChild_ != NULL && !treeCleaning_){
      candidateChild_ = NULL;
    }
    else
      CbcTree::pop();
  }
  /// Remove the best node from the heap and return it
  CbcNode *
  CbcProbedDiver::bestNode(double cutoff)
  {
#ifdef DIVE_DEBUG
    std::cout<<"CbcProbedDiver::bestNode"<<std::endl;
#endif
    if (nextOnBranch_ != NULL && !treeCleaning_) {
      if (nextOnBranch_->objectiveValue() < cutoff) {
        if (stop_diving_on_cutoff_ && nextOnBranch_->guessedObjectiveValue() >= cutoff) {
#ifdef DIVE_DEBUG
          printf("Stopping dive %g %g\n", nextOnBranch_->objectiveValue(), nextOnBranch_->guessedObjectiveValue());
#endif
          CbcTree::push(nextOnBranch_);
          nextOnBranch_ = NULL;
          //Also have to cleanup candidateChild_
          CbcTree::push(candidateChild_);
          candidateChild_ = NULL;
          return CbcTree::bestNode(cutoff);
        }
#ifdef DIVE_DEBUG
        printf("Diving on %i. obj=%g guesses=%g\n", nextOnBranch_, 
                nextOnBranch_->objectiveValue(), 
                nextOnBranch_->guessedObjectiveValue());
#endif
        CbcNode * ret_val = nextOnBranch_;
        nextOnBranch_ = NULL;
        return ret_val;
      }
      assert(true && "Should not get here.");
      CbcTree::push(nextOnBranch_);//For safety
      nextOnBranch_ = NULL;
    }
    else if(candidateChild_ != NULL && ! treeCleaning_){
#ifdef DIVE_DEBUG
      printf("brother was infeasible!!\n``");
#endif
      if(candidateChild_->objectiveValue() < cutoff) {
        if(stop_diving_on_cutoff_ && candidateChild_->guessedObjectiveValue() >= cutoff) {
#ifdef DIVE_DEBUG
          printf("Stopping dive %g %g\n", candidateChild_->objectiveValue(), candidateChild_->guessedObjectiveValue());
#endif
          CbcTree::push(candidateChild_);
          candidateChild_ = NULL;
          return CbcTree::bestNode(cutoff);
        }
#ifdef DIVE_DEBUG
        printf("Diving on %i. obj=%g guessed=%g\n", 
               candidateChild_,
               candidateChild_->objectiveValue(), 
               candidateChild_->guessedObjectiveValue());
#endif
        CbcNode * ret_val =  candidateChild_;
        candidateChild_ = NULL;
        return ret_val;
      }
      assert(true && "Should not get here.");
   }
   CbcNode * ret_val = CbcTree::bestNode(cutoff);
#ifdef DIVE_DEBUG
        printf("Picking top of the heap node %i", ret_val);
#endif
    return ret_val;
  }


  /** Test if empty. */
  bool CbcProbedDiver::empty()
  {
    return (CbcTree::empty() && (nextOnBranch_ == NULL) && (candidateChild_ == NULL) );
  }

  /*! \brief Prune the tree using an objective function cutoff
    if nextOnBranch_ exists we push it on the heap and call CbcTree function
  */
  void
  CbcProbedDiver::cleanTree(CbcModel * model, double cutoff, double & bestPossibleObjective)
  {
    if (nextOnBranch_ != NULL)
      CbcTree::push(nextOnBranch_);
    if (candidateChild_ != NULL)
      CbcTree::push(candidateChild_);
    nextOnBranch_=NULL;
    candidateChild_ = NULL;
    treeCleaning_ = true;
    CbcTree::cleanTree(model,cutoff, bestPossibleObjective);
    treeCleaning_ = false;
  }

  /// Get best possible objective function in the tree
  double
  CbcProbedDiver::getBestPossibleObjective()
  {
    double bestPossibleObjective = (nextOnBranch_ != NULL) ? nextOnBranch_->objectiveValue() : 1e100;
    if(candidateChild_ != NULL && candidateChild_->objectiveValue() < bestPossibleObjective )
       bestPossibleObjective = candidateChild_->objectiveValue();
    for (unsigned int i = 0 ; i < nodes_.size() ; i++) {
      if (nodes_[i] == NULL) continue;
      const double & obj = nodes_[i]->objectiveValue();
      if (obj < bestPossibleObjective) {
        bestPossibleObjective = obj;
      }
    }
    return bestPossibleObjective;
  }

  /// Initialize the method (get options)
  void
  CbcProbedDiver::initialize(BabSetupBase &b)
  {
    b.options()->GetBoolValue("stop_diving_on_cutoff", stop_diving_on_cutoff_,
        b.prefix());
  }

  /************************************************************************/
  /*                CbcDfsDiver methods                                    */
  /************************************************************************/
  /// Default constructor.
  CbcDfsDiver::CbcDfsDiver(): CbcTree(),
      treeCleaning_(false),
      dive_(),
      diveListSize_(0),
      divingBoardDepth_(-1),
      cutoff_(1e100),
      nBacktracks_(0),
      maxDepthBFS_(4),
      maxDiveBacktracks_(2),
      maxDiveDepth_(COIN_INT_MAX),
      mode_(Enlarge)
  {}

  ///Copy constructor.
  CbcDfsDiver::CbcDfsDiver(const CbcDfsDiver &rhs):CbcTree(rhs),
      treeCleaning_(rhs.treeCleaning_),
      dive_(rhs.dive_),
      diveListSize_(rhs.diveListSize_),
      divingBoardDepth_(rhs.divingBoardDepth_),
      cutoff_(rhs.cutoff_),
      nBacktracks_(rhs.nBacktracks_),
      maxDepthBFS_(rhs.maxDepthBFS_),
      maxDiveBacktracks_(rhs.maxDiveBacktracks_),
      maxDiveDepth_(rhs.maxDiveDepth_),
      mode_(rhs.mode_)
  {}
  /// Assignment operator.
  CbcDfsDiver &
  CbcDfsDiver::operator=(const CbcDfsDiver &rhs)
  {
    if (this != &rhs) {
      CbcTree::operator=(rhs);
      treeCleaning_ = rhs.treeCleaning_;
      dive_ = rhs.dive_;
      diveListSize_ = rhs.diveListSize_;
      divingBoardDepth_ = rhs.divingBoardDepth_;
      cutoff_ = rhs.cutoff_;
      nBacktracks_ = rhs.nBacktracks_;
      maxDepthBFS_ = rhs.maxDepthBFS_;
      maxDiveBacktracks_ = rhs.maxDiveBacktracks_;
      maxDiveDepth_ = maxDiveDepth_;
      mode_ = rhs.mode_;
    }
    return *this;
  }

  /// Destructor.
  CbcDfsDiver::~CbcDfsDiver()
  {}

  ///copy constructor.
  CbcTree *
  CbcDfsDiver::clone() const
  {
    return new CbcDfsDiver(*this);
  }

  ///Return top node (next node to process.*/
  CbcNode *
  CbcDfsDiver::top() const
  {
    if (treeCleaning_) return CbcTree::top();
#ifdef DIVE_DEBUG
    std::cout<<"CbcDfsDiver::top"<<std::endl;
#endif
    if (mode_ != CbcDfsDiver::FindSolutions) {
      assert(dive_.empty());
      CbcTree::top();
    }
    if (diveListSize_) {
      return dive_.front();
    }
    else return CbcTree::top();
  }

  /// Add node to the heap.
  void
  CbcDfsDiver::push(CbcNode * x)
  {
    if (treeCleaning_) return CbcTree::push(x);
#ifdef DIVE_DEBUG
    std::cout<<"CbcDfsDiver::push"<<std::endl;
#endif
    if (mode_ > CbcDfsDiver::FindSolutions) {
      CbcTree::push(x);
      assert(dive_.empty());
      return;
    }
    //Always push on dive;
    dive_.push_front(x);
    diveListSize_++;
#ifdef DIVE_DEBUG
    printf("diveListSize_ = %i == %u = dive_.size()\n",diveListSize_, dive_.size());
    assert(diveListSize_ == dive_.size());
#endif
  }

  /// Remove the top node of the heap.
  void
  CbcDfsDiver::pop()
  {
    if (treeCleaning_) return CbcTree::pop();
#ifdef DIVE_DEBUG
    std::cout<<"CbcDfsDiver::pop"<<std::endl;
#endif
    if (mode_ > CbcDfsDiver::FindSolutions) {
      assert(dive_.empty());
    }
    if (!dive_.empty()) {
      dive_.pop_front();
      diveListSize_--;
    }
    else
      CbcTree::pop();
  }

  /// Remove the best node from the heap and return it
  CbcNode *
  CbcDfsDiver::bestNode(double cutoff)
  {
#ifdef DIVE_DEBUG
    std::cout<<"CbcDfsDiver::bestNode"<<std::endl;
#endif
    if (treeCleaning_) return CbcTree::bestNode(cutoff);
#ifdef DIVE_DEBUG
    for (unsigned int i = 0 ; i < nodes_.size() ; i++) {
      if (nodes_[i]->objectiveValue() >= cutoff)
        std::cerr<<"CbcDfsDiver::bestNode"<<std::endl
        <<nodes_[i]->objectiveValue()<<", "<<cutoff<<std::endl;
      assert(nodes_[i]->objectiveValue() < cutoff);
    }
#endif
    if (mode_ == CbcDfsDiver::Enlarge) {
      if (diveListSize_ == 0)
        mode_ = CbcDfsDiver::FindSolutions;
      else {
        CbcNode * node = dive_.back();
        assert(node != NULL);
        if (node->depth() > maxDepthBFS_) {
          //Switch mode to Diving
          setComparisonMode(FindSolutions);
        }
        else {
          //pop and return node;
          dive_.pop_back();
          diveListSize_ --;
          return node;
        }
      }
    }
    if (mode_ != CbcDfsDiver::FindSolutions) {
      assert(dive_.empty());
      CbcTree::bestNode(cutoff);
    }
    assert(nBacktracks_ < maxDiveBacktracks_);
    CbcNode * node = NULL;
    while (diveListSize_ > 0) {
#ifdef DIVE_DEBUG
      std::cerr<<"CbcDfsDiver::bestNode"
      <<", examining node"<<std::endl;
#endif
      assert(!dive_.empty());
      node = dive_.front();
      dive_.pop_front();
      diveListSize_ --;
      assert(node);
      assert((node->depth() - divingBoardDepth_) <= maxDiveDepth_);
      if (node->objectiveValue() > cutoff) {//throw away node for now just put it on the heap as deleting a node is
        //more complicated than that (has to delete nodeInfo, cuts...)
#ifdef DIVE_DEBUG
        std::cout<<"CbcDfsDiver::bestNode"
        <<", node above cutoff"<<std::endl;
#endif
        CbcTree::push(node);
        node = NULL;
        nBacktracks_++;
      }
      else if (0 && node->guessedObjectiveValue() > cutoff) {//Put it on the real heap
#ifdef DIVE_DEBUG
        std::cout<<"CbcDfsDiver::bestNode"
        <<", node estimates "<<node->guessedObjectiveValue()<<"above cutoff"
        <<cutoff<<std::endl;
#endif
        CbcTree::push(node);
        nBacktracks_++;
        node = NULL;
      }
      else if ((node->depth() - divingBoardDepth_) > maxDiveDepth_) {//Put it on the real heap
#ifdef DIVE_DEBUG
        std::cout<<"CbcDfsDiver::bestNode"
        <<", node too deep"<<std::endl;
#endif
        CbcTree::push(node);
        nBacktracks_++;
        node = NULL;
      }
      else if (node->branchingObject()->numberBranchesLeft() < node->branchingObject()->numberBranches()) {//Backtracking
        nBacktracks_++;
#ifdef DIVE_DEBUG
        std::cout<<"CbcDfsDiver::bestNode"
        <<", backtracking"<<std::endl;
#endif
      }
      if (nBacktracks_ >= maxDiveBacktracks_) {//Push all the node in dive_ onto nodes_
#ifdef DIVE_DEBUG
        std::cout<<"CbcDfsDiver::bestNode"
        <<", maximum number of backtracks attained emptying dive_"<<std::endl;
#endif
        pushDiveOntoHeap(-COIN_DBL_MAX);
        if (node != NULL) CbcTree::push(node);
        node = NULL;
      }
      if (node != NULL)
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


  void CbcDfsDiver::pushDiveOntoHeap(double cutoff)
  {
    while (!dive_.empty() ){//&& dive_.front()->objectiveValue() >= cutoff) {
      CbcTree::push(dive_.front());
      dive_.pop_front();
      diveListSize_--;
    }
    for (std::list<CbcNode *>::iterator i = dive_.begin() ; i != dive_.end();
        i++) {
      assert(*i != NULL);
    }
  }
  /** Test if empty. */
  bool CbcDfsDiver::empty()
  {
    return (CbcTree::empty() && dive_.empty());
  }

  /*! \brief Prune the tree using an objective function cutoff
    if nextOnBranch_ exists we push it on the heap and call CbcTree function
  */
  void
  CbcDfsDiver::cleanTree(CbcModel * model, double cutoff, double & bestPossibleObjective)
  {
#ifdef DIVE_DEBUG
    std::cout<<"CbcDfsDiver::cleanTree"<<std::endl;
    std::cout<<"cutoff: "<<cutoff<<std::endl;
#endif
    pushDiveOntoHeap(cutoff);
    treeCleaning_ = true;
    CbcTree::cleanTree(model,cutoff, bestPossibleObjective);
    treeCleaning_ = false;
  }

  /// Get best possible objective function in the tree
  double
  CbcDfsDiver::getBestPossibleObjective()
  {
#ifdef DIVE_DEBUG
    std::cout<<"CbcDfsDiver::getBestPossibleObjective"<<std::endl;
#endif
    double bestPossibleObjective = CbcTree::empty() ? COIN_DBL_MAX : CbcTree::getBestPossibleObjective();
    for (std::list<CbcNode *>::iterator i = dive_.begin() ; i != dive_.end() ; i++) {
      if (*i == NULL) continue;
      const double & obj = (*i)->objectiveValue();
      if (obj < bestPossibleObjective) {
        bestPossibleObjective = obj;
      }
    }
    return bestPossibleObjective;
  }
//#ifdef COIN_HAS_BONMIN
  ///Register the options of the method.
  void
  CbcDfsDiver::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("Diving options", RegisteredOptions::UndocumentedCategory);
    roptions->AddLowerBoundedIntegerOption("max_backtracks_in_dive",
        "Set the number of backtracks in a dive when using dfs-dive tree search"
        "strategy.",
        0,5,
        "");
    roptions->setOptionExtraInfo("max_backtracks_in_dive",27);

    roptions->AddLowerBoundedIntegerOption("max_dive_depth",
        "When using dfs-dive search. Maximum depth to go to from the diving "
        "board (node where the diving started.",
        0,COIN_INT_MAX,
        "");
    roptions->setOptionExtraInfo("max_dive_depth",27);

  }

  /// Initialize the method (get options)
  void
  CbcDfsDiver::initialize(BabSetupBase &b)
  {
    b.options()->GetIntegerValue("max_dive_depth", maxDiveDepth_,b.prefix());
    b.options()->GetIntegerValue("max_backtracks_in_dive", maxDiveBacktracks_,b.prefix());
  }
//#endif

  /** Changes the mode of comparison of the tree for "safety reasons" if the mode really changes we always
      finish the current dive and put all the node back onto the heap.*/
  void
  CbcDfsDiver::setComparisonMode(ComparisonModes newMode)
  {
    if (newMode != mode_) {
      mode_ = newMode;
      //Empty heap
      pushDiveOntoHeap(-COIN_DBL_MAX);
      nBacktracks_ = maxDiveBacktracks_ -1;//Force to start a new dive
#ifdef DIVE_DEBUG
      std::cout<<"CbcDfsDiver::setComparisonMode"
      <<std::endl;
      switch (mode_) {
      case Enlarge:
        std::cout<<"Enlarge"<<std::endl;
        break;
      case FindSolutions:
        std::cout<<"FindSolutions"<<std::endl;
        break;
      case CloseBound:
        std::cout<<"CloseBound"<<std::endl;
        break;
      case LimitTreeSize:
        std::cout<<"LimitTreeSize"<<std::endl;
        break;
      }
#endif
      CbcTree::setComparison(*comparison_.test_); 
    }
  }



  // This allows any method to change behavior as it is called
  // after each solution
  bool
  DiverCompare::newSolution(CbcModel * model)
  {
    assert(diver_);
#ifdef DIVE_DEBUG
    std::cout<<"CbcDfsDiver::newSolution"
    <<std::endl;
    std::cout<<"Found "<<model->getSolutionCount()<<" solutions"<<std::endl;
#endif
    bool r_value = false;
    if (diver_->getComparisonMode() == CbcDfsDiver::Enlarge){
      diver_->setComparisonMode(CbcDfsDiver::FindSolutions);
      r_value = true;}
    if (model->getSolutionCount() >= numberSolToStopDive_ && diver_->getComparisonMode() == CbcDfsDiver::FindSolutions) {
      diver_->setComparisonMode(CbcDfsDiver::CloseBound);
      r_value = true;
    }
    return r_value;
  }

  /// This is test function
  bool
  DiverCompare::test (CbcNode * x, CbcNode * y)
  {
    assert(diver_);
    assert(comparisonDive_);
    assert(comparisonBound_);
    CbcDfsDiver::ComparisonModes mode = diver_->getComparisonMode();
    if (mode == CbcDfsDiver::FindSolutions) {
      return comparisonDive_->test(x,y);
    }
    else if (mode == CbcDfsDiver::CloseBound) {
      return comparisonBound_->test(x,y);
    }
    else if (mode == CbcDfsDiver::LimitTreeSize) {
      return comparisonDepth_.test(x,y);
    }
    else {
      throw CoinError("DiverCompare","test"," Unknown mode for comparison.");
    }
  }


  // This Also allows any method to change behavior as it is called
  // after each solution
  bool
  DiverCompare::newSolution(CbcModel * model,
      double objectiveAtContinuous,
      int numberInfeasibilitiesAtContinuous)
  { return false; }

  // This allows any method to change behavior as it is called
  // after every 1000 nodes.
  // Return true if want tree re-sorted
  bool
  DiverCompare::every1000Nodes(CbcModel * model,int numberNodes)
  {
    assert(diver_);
    if (numberNodes > numberNodesToLimitTreeSize_  && diver_->getComparisonMode() != CbcDfsDiver::LimitTreeSize) {
      diver_->setComparisonMode(CbcDfsDiver::LimitTreeSize);
      return true;
    }
    return false;
  }

}/* Ends namespace Bonmin.*/


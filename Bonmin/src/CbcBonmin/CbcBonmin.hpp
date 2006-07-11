#ifndef BonminBB_hpp
#define BonminBB_hpp

#include "BonminCbcParam.hpp"

class IpoptInterface;
class CbcObject;
/** Class which performs optimization of an MINLP stored in an IpoptInterface. */
class BonminBB
{
public:
  /** Integer optimization return codes.*/
  enum MipStatuses {FeasibleOptimal /** Optimum solution has been found and its optimality proved.*/,
      ProvenInfeasible /** Problem has been proven to be infeasible.*/,
      Feasible /** An integer solution to the problem has been found.*/,
      NoSolutionKnown/** No feasible solution to the problem is known*/};
  /** Constructor.*/
  BonminBB();
  /** destructor.*/
  virtual ~BonminBB();
  /** Perform a branch-and-bound on given IpoptInterface using passed parameters.*/
  void branchAndBound(IpoptInterface &ip,
      const BonminCbcParam&par);

  /**operator() performs the branchAndBound*/
  void operator()(IpoptInterface &ip, const BonminCbcParam& par)
  {
    branchAndBound(ip,par);
  }

  /** get the best solution known to the problem (is optimal if MipStatus is FeasibleOptimal).
      if no solution is known returns NULL.*/
  const double * bestSolution() const
  {
    return bestSolution_;
  }
  /** return objective value of the bestSolution */
  double bestObj() const
  {
    return bestObj_;
  }
  /** return Mip Status */
  MipStatuses mipStatus() const
  {
    return mipStatus_;
  }
  /** return the best known lower bound on the objective value*/
  double bestBound();
  /** return the total number of nodes explored.*/
  int numNodes() const
  {
    return numNodes_;
  }
  /** return the total number of iterations in the last mip solved.*/
  int iterationCount()
  {
    return mipIterationCount_;
  }
  /** returns the value of the continuous relaxation. */
  double continuousRelaxation()
  {
    return continuousRelaxation_;
  }

  /** virtual callback function to eventually modify objects for integer variable
      (replace with user set). This is called after CbcModel::findIntegers */
  virtual void replaceIntegers(CbcObject ** objects, int numberObjects)
  {};
private:
  /** Stores the solution of MIP. */
  double * bestSolution_;
  /** Status of the mixed integer program. */
  MipStatuses mipStatus_;
  /** objValue of MIP */
  double bestObj_;
  /** best known (lower) bound.*/
  double bestBound_;
  /** Continuous relaxation of the problem */
  double continuousRelaxation_;
  /** Number of nodes enumerated.*/
  int numNodes_;
  /** get total number of iterations in last mip solved.*/
  int mipIterationCount_;
};
#endif

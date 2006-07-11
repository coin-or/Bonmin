// (C) Copyright Carnegie Mellon University 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  05/26/2005

#ifndef IpCbcDummyHeuristic_HPP
#define IpCbcDummyHeuristic_HPP
#include "IpoptInterface.hpp"

#include "CbcHeuristic.hpp"

class  IpCbcDummyHeuristic : public CbcHeuristic
{
public:
  /// Default constructor
  IpCbcDummyHeuristic(IpoptInterface * si = NULL);
  /// Usefull constructor
  IpCbcDummyHeuristic(CbcModel &model, IpoptInterface * si = NULL);
  ///Copy constructor
  IpCbcDummyHeuristic( const IpCbcDummyHeuristic &copy):
      CbcHeuristic(copy),
      nlp_(copy.nlp_),
      knowsSolution(copy.knowsSolution)
  {}
  /// Assign an IpoptInterface
  void assignInterface(IpoptInterface * si);
  /// heuristic method
  virtual int solution(double &solutionValue, double *betterSolution);
  virtual int solution(double &solutionValue, double *betterSolution, OsiCuts & cs)
  {
    return solution(solutionValue, betterSolution);
  }
  virtual CbcHeuristic * clone()const
  {
    return new IpCbcDummyHeuristic(*this);
  }
  virtual void resetModel(CbcModel*)
  {}
private:
  /// Pointer to the Ipopt interface
  IpoptInterface * nlp_;
  /// Do I have a solution?
  bool knowsSolution;
};
#endif

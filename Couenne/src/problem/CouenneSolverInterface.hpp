/*
 * Name:    CouenneSolverInterface.hpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: OsiSolverInterface with a pointer to a CouenneProblem object
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNESOLVERINTERFACE_HPP
#define COUENNESOLVERINTERFACE_HPP

#include <CoinFinite.hpp>
#include <OsiClpSolverInterface.hpp>
#include <CouenneProblem.hpp>

/// TODO: 
/// 1) bound tightening before solver
/// 2) re-implement OSI::isInteger with problem_ -> whatever -> isInteger ()

class CouenneSolverInterface: public OsiClpSolverInterface {

private:

  CouenneProblem *problem_;

public:

  CouenneSolverInterface () {}
  ~CouenneSolverInterface () {}

  CouenneProblem Problem ()
  {return problem_;}
}

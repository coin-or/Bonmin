// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/13/2007

#ifndef BonBaseOptions_H
#define BonBaseOptions_H
namespace Bonmin {
namespace BaseOptions {
  /** Type of algorithms which can be used.*/
  enum Algorithm{
    B_BB=0/** Bonmin's Branch-and-bound.*/,
    B_OA=1/** Bonmin's Outer Approximation Decomposition.*/,
    B_QG=2/** Bonmin's Quesada & Grossmann branch-and-cut.*/,
    B_Hyb=3/** Bonmin's hybrid outer approximation.*/,
    B_Couenne=4/** Bonmin's and Couenne spatial branch-and-bound.*/
  };
  
  /** Solvers for solving nonlinear programs.*/
  enum Solver{
    Ipopt=0 /** <a href="http://projects.coin-or.org/Ipopt">
    Ipopt </a> interior point algorithm.*/,
    FilterSQP /** <a href="http://www-unix.mcs.anl.gov/~leyffer/solvers.html"> filterSQP </a> Sequential Quadratic Programming algorithm.*/
  };
}
}
#endif
/*
 * Name: CouenneCglCutGenerator.h
 * Author: Pietro Belotti
 * Purpose: define a class of pools of convexification OsiRowCuts
 *
 * (C) Pietro Belotti, all rights reserved. 
 * This file is licensed under the Common Public License.
 */

#ifndef COUENNE_CUT_GENERATOR_H
#define COUENNE_CUT_GENERATOR_H

#include <iostream>

#include <CglCutGenerator.hpp>
#include <OsiRowCut.hpp>
#include <CouenneProblem.h>

#include <OsiSolverInterface.hpp>
#include <OsiClpSolverInterface.hpp>


struct ASL_pfgh;


// Cut Generator for linear convexifications

class CouenneCutGenerator: public CglCutGenerator {

 protected:

  // number of cuts currently in pool
  int ncuts_;

  // array of pointers to linearization cuts
  OsiRowCut **pool_;

  // for use with Bonmin. It is NULL when used normally
  OsiCuts *bonCs_;

  // similarly, a fictitious object to make calls to generateCuts 
  OsiSolverInterface *bonOs_;

  // has generateCuts been called yet?
  mutable bool firstcall_;

  // should we add the violated cuts only (true), or all of them (false)?
  bool addviolated_;

  // what kind of sampling should be performed?
  enum conv_type convtype_;

  // how many cuts should be added for each function?
  int nSamples_;

  // pointer to symbolic representation of constraint, variables, and
  // bounds
  CouenneProblem *problem_;

 public:

  // constructor
  CouenneCutGenerator  (const struct ASL_pfgh * = NULL, 
			bool = false, 
			enum conv_type = UNIFORM_GRID, 
			int = 2);
  // copy constructor
  CouenneCutGenerator  (const CouenneCutGenerator &);

  // destructor
  ~CouenneCutGenerator ();

  // clone method (necessary for the abstract CglCutGenerator class)
  CouenneCutGenerator *clone () const;

  // return pointer to symbolic problem
  CouenneProblem *Problem () const
    {return problem_;}

  // get methods
  int         getncuts   ()      const {return ncuts_;}
  OsiRowCut  *getCut     (int i) const {return pool_ [i];}

  const CouNumber   X    (int i) const {return problem_ -> X  (i);}

  const CouNumber   &Lb  (int i) {return problem_ -> Lb (i);}
  const CouNumber   &Ub  (int i) {return problem_ -> Ub (i);}

  // get arrays
  const CouNumber   *X   ()      {return problem_ -> X  ();}
  const CouNumber   *Lb  ()      {return problem_ -> Lb ();}
  const CouNumber   *Ub  ()      {return problem_ -> Ub ();}

  int getnvars () const {return problem_ -> nVars () + 
			        problem_ -> nAuxs ();} 

  // Solver interface (used in copy constructor)
  OsiSolverInterface *getBonOs () const
    {return bonOs_;}

  // Cutset (used in copy constructor)
  OsiCuts *getBonCs () const
    {return bonCs_;}

  // has generateCuts been called yet?
  bool isFirst () const 
    {return firstcall_;}

  // should we add the violated cuts only (true), or all of them (false)?
  bool addViolated () const
    {return addviolated_;}

  // get convexification type (see CouenneTypes.h)
  enum conv_type ConvType () const
    {return convtype_;}

  // get number of convexification samples
  int nSamples () const
    {return nSamples_;}

  // the main CglCutGenerator
  void generateCuts (const OsiSolverInterface &, 
		     OsiCuts &, 
		     const CglTreeInfo = CglTreeInfo ()) const;

  // update linearization cut pool (parameters are current point,
  // current lower-, and current upper bound)
  int updateConv (CouNumber *, CouNumber *, CouNumber *);

  // create cut and check violation
  OsiRowCut *createCut (CouNumber, // rhs
			int,       // sign: -1: <=, 0: =, +1: >=
			int,    CouNumber,    // index, coeff of first term
			int=-1, CouNumber=0., //              of second term 
			                      // (-1 means don't consider)
			int=-1, CouNumber=0.) const; //       of third term

  // add general linear envelope to convex function, given its
  // variables' indices, the (univariate) function and its first
  // derivative
  void addEnvelope (OsiCuts &,
		    int,
		    unary_function, unary_function, 
		    int, int, 
		    CouNumber, CouNumber, CouNumber) const;

  // add half-plane through (x1,y1) and (x2,y2) -- resp. 4th, 5th,
  // 6th, and 7th argument
  void addSegment (OsiCuts &, int, int,
		   CouNumber, CouNumber, 
		   CouNumber, CouNumber, int) const;

  // add half-plane corresponding to tangent in given point
  void addTangent (OsiCuts &, int, int, CouNumber, CouNumber, CouNumber, int) const;
};

#endif

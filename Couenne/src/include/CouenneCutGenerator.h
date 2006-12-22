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
#include <CouenneProblem.h>


struct ASL_pfgh;


// Cut Generator for linear convexifications

class CouenneCutGenerator: public CglCutGenerator {

 protected:

  // number of cuts currently in pool
  int ncuts_;

  // array of pointers to linearization cuts
  OsiRowCut **pool_;

  // for compatibility with Bonmin. It is NULL when used normally
  //  OsiCuts *bonCs_;

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
  // destructor
  ~CouenneCutGenerator () {}

  // clone method (necessary for the abstract CglCutGenerator class)
  CouenneCutGenerator *clone () const 
    {return NULL;}

  // return pointer to symbolic problem
  CouenneProblem *Problem ()  
    {return problem_;}

  // get methods
  int         getncuts ()      const {return ncuts_;}
  OsiRowCut  *getCut   (int i) const {return pool_ [i];}

  const CouNumber   &X       (int i)  {return problem_ -> X  (i);}
  const CouNumber   &Lb      (int i)  {return problem_ -> Lb (i);}
  const CouNumber   &Ub      (int i)  {return problem_ -> Ub (i);}

  const CouNumber   *X       ()       {return problem_ -> X  ();}
  const CouNumber   *Lb      ()       {return problem_ -> Lb ();}
  const CouNumber   *Ub      ()       {return problem_ -> Ub ();}

  int getnvars () const {return problem_ -> nVars () + 
			        problem_ -> nAuxs ();} 

  // has generateCuts been called yet?
  bool isFirst () const 
    {return firstcall_;}

  // should we add the violated cuts only (true), or all of them (false)?
  bool addViolated () const
    {return addviolated_;}

  // get or set convexification type (see CouenneTypes.h)
  enum conv_type ConvType () const
    {return convtype_;}

  // get or set convexification type (see CouenneTypes.h)
  int nSamples () const
    {return nSamples_;}

  // the main CglCutGenerator
  void generateCuts (const OsiSolverInterface &, 
		     OsiCuts &, 
		     const CglTreeInfo) const;

  // for compatibility with Bonmin
  //  int generateCuts (const OsiSolverInterface &, 
  //		    OsiRowCut **&);

  // (re-)initializes the pool  
  void cleanup ();

  // update linearization cut pool (parameters are current point,
  // current lower-, and current upper bound)
  int updateConv (CouNumber *, CouNumber *, CouNumber *);

  // updates auxiliary variables' bounds
  void updateBounds (CouNumber *, CouNumber *);
};


// add half-plane corresponding to tangent in given point

void addTangent (OsiCuts &, int, int, CouNumber, CouNumber, CouNumber, int);


// add half-plane defined by two points (x1,y1) and (x2, y2), and
// sign. Sign is only valid if x1 < x2

inline void addSegment (OsiCuts &cs, int wi, int xi, 
			CouNumber x1, CouNumber y1, 
			CouNumber x2, CouNumber y2, int sign) { 

  if (x2-x1 > COUENNE_EPS) {
    printf ("Segment --> ");
    addTangent (cs, wi, xi, x1, y1, (y2-y1) / (x2-x1), sign);
  }
}

#endif

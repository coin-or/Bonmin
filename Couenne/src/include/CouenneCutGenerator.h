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
#include <BonOaDecBase.hpp>
#include <OsiRowCut.hpp>
#include <CouenneTypes.h>

#include <OsiSolverInterface.hpp>
#include <OsiClpSolverInterface.hpp>

class CouenneProblem;

struct ASL_pfgh;


/// Cut Generator for linear convexifications

class CouenneCutGenerator: public Bonmin::OaDecompositionBase {

 protected:

  /// has generateCuts been called yet?
  mutable bool firstcall_;

  /// should we add the violated cuts only (true), or all of them (false)?
  bool addviolated_;

  /// what kind of sampling should be performed?
  enum conv_type convtype_;

  /// how many cuts should be added for each function?
  int nSamples_;

  /// pointer to symbolic repr. of constraint, variables, and bounds
  CouenneProblem *problem_;

  /// number of cuts generated at the first call
  mutable int nrootcuts_;

  /// total number of cuts generated 
  mutable int ntotalcuts_;

  /// Record obj value at final point of CouenneConv.
  mutable double objValue_;

  /// nonlinear solver interface as used within Bonmin (used at first
  /// Couenne pass of each b&b node
  Bonmin::OsiTMINLPInterface *nlp_;

 public:

  /// constructor
  CouenneCutGenerator  (Bonmin::OsiTMINLPInterface * = NULL,
			const struct ASL_pfgh * = NULL, 
			bool = false, 
			enum conv_type = UNIFORM_GRID, 
			int = 2);
  /// copy constructor
  CouenneCutGenerator  (const CouenneCutGenerator &);

  /// destructor
  ~CouenneCutGenerator ();

  /// clone method (necessary for the abstract CglCutGenerator class)
  CouenneCutGenerator *clone () const
  {return new CouenneCutGenerator (*this);}

  /// return pointer to symbolic problem
  CouenneProblem *Problem () const
    {return problem_;}

  const CouNumber   X    (int);

  const CouNumber   &Lb  (int);
  const CouNumber   &Ub  (int);

  /// get arrays
  const CouNumber   *X   ();
  const CouNumber   *Lb  ();
  const CouNumber   *Ub  ();

  int getnvars () const;

  /// has generateCuts been called yet?
  bool isFirst () const 
    {return firstcall_;}

  /// should we add the violated cuts only (true), or all of them (false)?
  bool addViolated () const
    {return addviolated_;}

  /// get convexification type (see CouenneTypes.h)
  enum conv_type ConvType () const
    {return convtype_;}

  /// get number of convexification samples
  int nSamples () const
    {return nSamples_;}

  /// the main CglCutGenerator
  void generateCuts (const OsiSolverInterface &, 
		     OsiCuts &, 
		     const CglTreeInfo = CglTreeInfo ()) const;

  /// create cut and check violation
  OsiRowCut *createCut (CouNumber, // rhs
			int,       // sign: -1: <=, 0: =, +1: >=
        	 		   // index, coeff  (index = -1 means "don't consider") 
			int,    CouNumber,    // of first  term
			int=-1, CouNumber=0., // of second term 
			int=-1, CouNumber=0., // of third  term
			bool = false) const;  // is it a global cut? No, by default

  /// add general linear envelope to convex function, given its
  /// variables' indices, the (univariate) function and its first
  /// derivative
  void addEnvelope (OsiCuts &,
		    int,
		    unary_function, unary_function, 
		    int, int, 
		    CouNumber, CouNumber, CouNumber,
		    bool = false) const;

  /// add half-plane through (x1,y1) and (x2,y2) -- resp. 4th, 5th,
  /// 6th, and 7th argument
  void addSegment (OsiCuts &, int, int,
		   CouNumber, CouNumber, 
		   CouNumber, CouNumber, int) const;

  /// virtual method which performs the OA algorithm by modifying lp and nlp.
  virtual double performOa (OsiCuts & cs,           solverManip &nlpManip, 
			    solverManip &lpManip,   SubMipSolver *& subMip, 
			    OsiBabSolver * babInfo, double &cutoff) const
    {throw -1;}

  /// virtual method to decide if local search is performed
  virtual bool doLocalSearch () const {return 0;}

  /// bound tightening
  virtual int tightenBounds (const OsiSolverInterface &, char *) const;
};

#endif

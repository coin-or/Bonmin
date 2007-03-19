/*
 * Name: CouenneCglCutGenerator.cpp
 * Author: Pietro Belotti
 * Purpose: define a class of convexification procedures 
 *
 * (C) Pietro Belotti, all rights reserved. 
 * This file is licensed under the Common Public License.
 */

// for profiling
#include <sys/time.h>

#include <OsiRowCut.hpp>
#include <BonOaDecBase.hpp>
#include <CglCutGenerator.hpp>
#include <CoinHelperFunctions.hpp>

#include <CouennePrecisions.h>
#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>


/// constructor

CouenneCutGenerator::CouenneCutGenerator (Bonmin::OsiTMINLPInterface *nlp,
					  const struct ASL *asl, 
					  bool addviolated,
					  enum conv_type convtype, int nSamples):

  OaDecompositionBase (nlp, NULL, NULL, 0,0,0),

  firstcall_      (true),
  addviolated_    (addviolated),
  convtype_       (convtype),
  nSamples_       (nSamples),
  problem_        (NULL),
  nrootcuts_      (0),
  ntotalcuts_     (0),
  objValue_       (- DBL_MAX),
  nlp_            (nlp),
  BabPtr_         (NULL) {

  if (!asl) return;

  problem_ = new CouenneProblem;

  double now = CoinCpuTime ();

  problem_ -> readnl      (asl);

  if ((now = (CoinCpuTime () - now)) > 10.)
    printf ("reading time: %.3fs\n", now);

  now = CoinCpuTime ();
  //  problem_ -> print (std::cout);
  //  printf ("======================================\n");
  problem_ -> standardize ();

  if ((now = (CoinCpuTime () - now)) > 10.)
  printf ("standardization time: %.3fs\n", now);

  septime_ = now;

  //  problem_ -> print (std::cout);
  //  printf ("======================================\n");
  //  problem_ -> writeMod ("extended.mod");
}


/// destructor

CouenneCutGenerator::~CouenneCutGenerator ()
{delete problem_;}


/// copy constructor

CouenneCutGenerator::CouenneCutGenerator (const CouenneCutGenerator &src):

  OaDecompositionBase (src),

  firstcall_   (src. firstcall_),
  addviolated_ (src. addviolated_), 
  convtype_    (src. convtype_), 
  nSamples_    (src. nSamples_),
  problem_     (src. problem_ -> clone ()),
  nrootcuts_   (src. nrootcuts_),
  ntotalcuts_  (src. ntotalcuts_),
  objValue_    (src. objValue_),
  nlp_         (src. nlp_),
  BabPtr_      (src. BabPtr_),
  septime_     (src. septime_)
{}


/// add half-space through two points (x1,y1) and (x2,y2)

void CouenneCutGenerator::addSegment (OsiCuts &cs, int wi, int xi, 
				      CouNumber x1, CouNumber y1, 
				      CouNumber x2, CouNumber y2, int sign) const { 

  OsiRowCut *cut = NULL;

  if (fabs (x2-x1) < COUENNE_EPS)
    // well, we suppose that y1=y2...
    cut = createCut (y2, 0, wi, CouNumber (1.));
  else {

    CouNumber oppslope = (y1-y2) / (x2-x1);

    cut = createCut (y1 + oppslope * x1, sign, 
		     wi, CouNumber (1.), 
		     xi, oppslope);
  }

  if (cut) 
    cs.insert (cut);
}


/// methods to get values of variables and bounds
const CouNumber    CouenneCutGenerator::X   (int i) {return problem_ -> X  (i);}
const CouNumber   &CouenneCutGenerator::Lb  (int i) {return problem_ -> Lb (i);}
const CouNumber   &CouenneCutGenerator::Ub  (int i) {return problem_ -> Ub (i);}


/// get arrays
const CouNumber   *CouenneCutGenerator::X   ()      {return problem_ -> X  ();}
const CouNumber   *CouenneCutGenerator::Lb  ()      {return problem_ -> Lb ();}
const CouNumber   *CouenneCutGenerator::Ub  ()      {return problem_ -> Ub ();}

int                CouenneCutGenerator::getnvars () const {return problem_ -> nVars () + 
					                          problem_ -> nAuxs ();} 

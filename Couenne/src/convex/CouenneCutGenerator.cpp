/*
 * Name: CouenneCglCutGenerator.cpp
 * Author: Pietro Belotti
 * Purpose: define a class of convexification procedures 
 *
 * (C) Pietro Belotti, all rights reserved. 
 * This file is licensed under the Common Public License.
 */

// for profiling
//AW 2007-06-11: Taken out this since it doesn't compile with MSVC++
//#include <sys/time.h>

#include <OsiRowCut.hpp>
#include <BonOaDecBase.hpp>
#include <CglCutGenerator.hpp>

#include <CouennePrecisions.h>
#include <CouenneProblem.hpp>
#include <CouenneCutGenerator.hpp>


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
  septime_        (0),
  objValue_       (- DBL_MAX),
  nlp_            (nlp),
  BabPtr_         (NULL),
  infeasNode_     (false) {

  problem_ = new CouenneProblem (asl);
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
  septime_     (src. septime_),
  objValue_    (src. objValue_),
  nlp_         (src. nlp_),
  BabPtr_      (src. BabPtr_),
  infeasNode_  (src. infeasNode_)
{}


/// add half-space through two points (x1,y1) and (x2,y2)

#define MAX_SLOPE 1e3

int CouenneCutGenerator::addSegment (OsiCuts &cs, int wi, int xi, 
				     CouNumber x1, CouNumber y1, 
				     CouNumber x2, CouNumber y2, int sign) const { 

  if (fabs (x2-x1) < COUENNE_EPS) {
    if (fabs (y2-y1) > MAX_SLOPE * COUENNE_EPS)
      printf ("warning, discontinuity of %e over an interval of %e\n", y2-y1, x2-x1);
    else return createCut (cs, y2, 0, wi, 1.);
  }

  //CouNumber oppslope = (y1-y2) / (x2-x1);

  CouNumber dx = x2-x1, dy = y2-y1;

  //  return createCut (cs, y1 + oppslope * x1, sign, wi, 1., xi, oppslope);
  return createCut (cs, y1*dx - dy*x1, (dx>0) ? sign : -sign, wi, dx, xi, -dy);
}


/// add tangent at (x,w) with given slope

int CouenneCutGenerator::addTangent (OsiCuts &cs, int wi, int xi, 
				     CouNumber x, CouNumber w, 
				     CouNumber slope, int sign) const { 

  return createCut (cs, w - slope * x, sign, wi, 1., xi, - slope);
}


/// total number of variables (original + auxiliary) of the problem
const int CouenneCutGenerator::getnvars () const
  {return problem_ -> nVars ();} 

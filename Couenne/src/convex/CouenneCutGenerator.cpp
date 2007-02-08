/*
 * Name: CouenneCglCutGenerator.cpp
 * Author: Pietro Belotti
 * Purpose: define a class of convexification procedures 
 *
 * (C) Pietro Belotti, all rights reserved. 
 * This file is licensed under the Common Public License.
 */

#include <OsiRowCut.hpp>
#include <CglCutGenerator.hpp>
#include <CoinHelperFunctions.hpp>

#include <CouenneCutGenerator.h>


// constructor

CouenneCutGenerator::CouenneCutGenerator (const ASL_pfgh *asl, bool addviolated,
					  enum conv_type convtype, int nSamples):
  CglCutGenerator (),
  ncuts_          (0),
  pool_           (NULL),
  bonCs_          (NULL),
  bonOs_          (NULL),
  firstcall_      (true),
  addviolated_    (addviolated),
  convtype_       (convtype),
  nSamples_       (nSamples),
  problem_        (NULL) {

  if (!asl) return;

  problem_ = new CouenneProblem;

  problem_ -> readnl      (asl);
  //  problem_ -> print (std::cout);
  //  printf ("======================================\n");
  problem_ -> standardize ();
  //  problem_ -> print (std::cout);
}


// destructor

CouenneCutGenerator::~CouenneCutGenerator () {

  if (bonCs_) {
    delete bonCs_;
    delete bonOs_;
  }

  if (!pool_) 
    free (pool_);

  delete problem_;
}


// clone method

CouenneCutGenerator *CouenneCutGenerator::clone () const
  {return new CouenneCutGenerator (*this);}


// copy constructor

CouenneCutGenerator::CouenneCutGenerator (const CouenneCutGenerator &src):
  CglCutGenerator (),
  ncuts_       (src.getncuts ()),
  pool_        (new OsiRowCut * [ncuts_]),
  bonCs_       (new OsiCuts (*(src.getBonCs ()))),
  bonOs_       (src.getBonOs () -> clone ()),
  firstcall_   (src.isFirst ()),
  addviolated_ (src. addViolated ()), 
  convtype_    (src. ConvType ()), 
  nSamples_    (src. nSamples ()),
  problem_     (src.Problem() -> clone ()) {

  for (int i=0; i<ncuts_; i++)
    pool_ [i] = src. getCut (i) -> clone ();
}


// add half-space through two points (x1,y1) and (x2,y2)

void CouenneCutGenerator::addSegment (OsiCuts &cs, int wi, int xi, 
				      CouNumber x1, CouNumber y1, 
				      CouNumber x2, CouNumber y2, int sign) const { 

  if (fabs (x2-x1) < COUENNE_EPS) {
    fprintf (stderr, "CouenneCutGenerator::addSegment(): warning, x1=x2\n");
    return;
  }

  CouNumber oppslope = (y1-y2) / (x2-x1);

  OsiRowCut *cut = createCut (y1 + oppslope * x1, sign, 
			      wi, CouNumber (1.), 
			      xi, oppslope);

  //  if (cut) {printf ("Segment: "); cut -> print ();}

  if (cut) 
    cs.insert (cut);
}

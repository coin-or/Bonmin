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
  problem_ -> standardize ();     
}


// destructor

CouenneCutGenerator::~CouenneCutGenerator () {

  if (bonCs_) {
    delete bonCs_;
    delete bonOs_;
  }
}


// clone method

CouenneCutGenerator *CouenneCutGenerator::clone () const
  {return new CouenneCutGenerator (*this);}

/*
  CouenneCutGenerator *dolly = new CouenneCutGenerator 
    (NULL, addviolated_, convtype_, nSamples_);

  dolly -> setIsFirst (firstcall_);
  dolly -> setBonCs   (bonCs_);
  dolly -> setBonOs   (bonOs_);

  dolly -> setProblem (problem_ -> clone ());
  dolly -> copyPool   (ncuts_, pool_);

  return dolly;
}
*/

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


// (re-)initializes cut pool

void CouenneCutGenerator::cleanup () {

  while (ncuts_--)
    delete pool_ [ncuts_];
  free (pool_);
  pool_ = NULL;
}


// Another version, to be called from Bonmin
/*
int CouenneCutGenerator::generateCuts (const OsiSolverInterface &si, OsiRowCut **&cuts) {

  int ncuts = 0;

  if (bonCs_ == NULL)
    bonCs_ = new OsiCuts;
  else {
    ncuts = bonCs_ -> sizeRowCuts ();
    for (int i=0; i<ncuts; i++)
      bonCs_ -> eraseRowCut (i);
  }

  generateCuts (si, *bonCs_);

  ncuts = bonCs_ -> sizeRowCuts ();

  cuts = new OsiRowCut * [ncuts];

  for (int i=0; i<ncuts; i++)
    cuts [i]= bonCs_ -> rowCutPtr (i);

  if (firstcall_)
    firstcall_ = false;

  return ncuts;
}
*/

void CouenneCutGenerator::addSegment (OsiCuts &cs, int wi, int xi, 
		 CouNumber x1, CouNumber y1, 
		 CouNumber x2, CouNumber y2, int sign) const { 

  CouNumber oppslope = (y1-y2)/(x2-x1);

  OsiRowCut *cut = createCut (y1 + oppslope * x1, sign, 
			      wi, CouNumber (1.), 
			      xi, oppslope);

  //  if (cut) {printf ("Segment: "); cut -> print ();}

  if (cut) 
    cs.insert (cut);
}

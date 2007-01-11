/*
 * Name: CouenneCglCutGenerator.cpp
 * Author: Pietro Belotti
 * Purpose: define a class of convexification procedures 
 *
 * (C) Pietro Belotti, all rights reserved. 
 * This file is licensed under the Common Public License.
 */

#include <CglCutGenerator.hpp>
#include <CouenneCutGenerator.h>


// constructor

CouenneCutGenerator::CouenneCutGenerator (const ASL_pfgh *asl, bool addviolated,
					  enum conv_type convtype, int nSamples):
  CglCutGenerator (),
  ncuts_          (0),
  pool_           (NULL),
  bonCs_          (NULL),
  firstcall_      (true),
  addviolated_    (addviolated),
  convtype_       (convtype),
  nSamples_       (nSamples),
  problem_        (new CouenneProblem) {

  if (!asl) return;

  problem_ -> readnl      (asl);  

  printf ("-- Original\n");     
  problem_ -> print (std::cout); 

  problem_ -> standardize ();     

  printf ("--------------------------\n");
  //  printf ("-- Standardized\n"); problem_ -> print (std::cout);
  //  problem_ -> convexify   ();     
  //  printf ("-- Convexified\n");  problem_ -> print (std::cout);
}


// destructor

CouenneCutGenerator::~CouenneCutGenerator () {

  if (bonCs_) {
    delete bonCs_;
    delete bonOs_;
  }
}


// (re-)initializes cut pool

void CouenneCutGenerator::cleanup () {

  while (ncuts_--)
    delete pool_ [ncuts_];
  free (pool_);
  pool_ = NULL;
}


// a convexifier cut generator

void CouenneCutGenerator::generateCuts (const OsiSolverInterface &si, 
					OsiCuts &cs, 
					const CglTreeInfo info) const {
  int ncols = si.getNumCols ();

  CouNumber *x, *l, *u;

  if (firstcall_) {

    // OsiSolverInterface is empty yet, no information can be obtained
    // on variables or bounds -- and none is needed since our
    // constructor populated *problem_ with variables and bounds. We
    // only need to 
  }
  else {

    // Retrieve, from si, variable and bounds of all variables, if not
    // firstcall, otherwise only those of the original ones Update
    // expression structure with x, l, u

    const CouNumber *xc = si.getColSolution (), 
      *lc = si.getColLower (),
      *uc = si.getColUpper ();

    x = new CouNumber [ncols];
    l = new CouNumber [ncols];
    u = new CouNumber [ncols];

    for (int i=ncols; i--;) {

      x [i] = xc [i];
      l [i] = lc [i];
      u [i] = uc [i];
    }

    printf ("!!!!-------x\n");
    for (int i=0; i < getnvars (); i++)
      printf ("%3d %12.9f\n", i, x [i]);

    problem_ ->  update (x,l,u);
    expression:: update (x,l,u);
  }

  // For each auxiliary variable replacing the original constraints,
  // check if corresponding bounds are violated, and add cut to cs

  if (firstcall_) {

    int nnlc = problem_ -> nNLCons ();

    for (int i=0; i<nnlc; i++) {

      CouenneConstraint *con = problem_ -> NLCon (i);

      // for constraint lb <= w <= ub, compute actual values of lb, w,
      // and ub

      //      CouNumber body = con -> Body () -> Value ();
      CouNumber lb   = con -> Lb   () -> Value ();
      CouNumber ub   = con -> Ub   () -> Value ();

      // if there exists violation, add constraint

      OsiRowCut *orc   = new OsiRowCut;
      CouNumber *coeff = new CouNumber [1];
      int       *index = new int       [1];

      coeff [0] = 1;
      index [0] = con -> Body () -> Index ();

      orc -> setRow (1, index, coeff);

      if (lb > - COUENNE_INFINITY + 1) orc -> setLb  (lb);
      if (ub <   COUENNE_INFINITY - 1) orc -> setUb  (ub);


      printf ("con %d: ", i);
      if (lb > - COUENNE_INFINITY) 
	printf ("%.4f <= ", lb);
      printf ("w_%d", index [0]);
      if (ub <   COUENNE_INFINITY) 
	printf (" <= %.4f\n", ub);

      printf ("1st... cut: "); orc -> print ();

      cs.insert (orc);
    }
  }

  // For each auxiliary variable, create cut (or set of cuts) violated
  // by current point and add it to cs

  for (int i = 0; i<problem_ -> nAuxs (); i++)
    problem_ -> Aux (i) -> generateCuts (si, cs, this);

  // end of generateCuts

  if (firstcall_) 
    firstcall_ = false;
  else {

    delete [] x;
    delete [] l;
    delete [] u;
  }
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

  if (cut) {printf ("Segment: "); cut -> print ();}

  if (cut) 
    cs.insert (cut);
}

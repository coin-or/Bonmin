/*
 * Name:    exprAux.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the class of auxiliary variables
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneCutGenerator.h>
#include <CouenneTypes.h>
#include <expression.h>
#include <exprAux.h>
#include <exprOp.h>
#include <exprUnary.h>
#include <exprVar.h>
#include <exprBound.h>

#include <CouenneProblem.h>

// auxiliary expression Constructor

exprAux::exprAux (expression *image, int index, int rank): 

  exprVar  (index),
  image_   (image),
  rank_    (rank),
  multiplicity_ (1) {

  // do this later, in standardize()
  //  image_ -> getBounds (lb_, ub_);
  getBounds (lb_, ub_);
}


// generate cuts for expression associated with this auxiliary

void exprAux::generateCuts (const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  static bool first_draw = true;
  static CouNumber maxY = -COUENNE_INFINITY,
                   minY =  COUENNE_INFINITY;

  int j = cs.sizeRowCuts ();
  CouNumber l;

  if ((!(cg -> isFirst ())) && 
      (fabs ((l = expression::Ubound (varIndex_)) - 
                  expression::Lbound (varIndex_)) < COUENNE_EPS))
    cg -> createCut (cs, l, 0, varIndex_, 1.);
  else image_ -> generateCuts (this, si, cs, cg);

  //  if (!(cg -> isFirst ())) 
  //  if (j < cs.sizeRowCuts ())
  if (0)
    {
      printf ("----------------Generated cut for "); 
      print (std::cout);  printf (" := ");
      image_ -> print (std::cout); 
      printf("\n");
      for (int jj=j; jj < cs.sizeRowCuts ();jj++)
	cs.rowCutPtr (jj) -> print ();
    }

  //////////////////////////////////////////////////////////////

  if (0) { // [cool!] print graph-readable output for displaying
           // inequalities on a Cartesian plane

    if (1 || (image_ -> code () == COU_EXPRSIN) || 
	(image_ -> code () == COU_EXPRPOW) || 
	(image_ -> code () == COU_EXPREXP) || 
	(image_ -> code () == COU_EXPRLOG) || 
	(image_ -> code () == COU_EXPRCOS)) {

      printf (" ==> "); print (std::cout); printf ("\n");

      expression *lbe, *ube;

      expression *indep = image_ -> Argument ();

      if (!indep) 
	indep = image_ -> getFixVar ();

      int xi = indep -> Index ();
      printf ("looking into w_%d = f (x_%d)\n", varIndex_, xi);

      indep -> getBounds (lbe, ube);

      CouNumber lb = (*lbe) (),
	        ub = (*ube) ();

      delete lbe;
      delete ube;

      if (xi >= 0) {

	CouNumber curx = expression::Variable (xi);

#define N_STEPS 100

	// plot function

	if (first_draw) {

	  first_draw = false;

	  for (CouNumber x = lb; 
	       x <= ub + COUENNE_EPS; 
	       x += ((ub - lb) / N_STEPS)) {

	    cg -> Problem () -> X () [xi] = x;
	  
	    CouNumber y = (*image_) ();

	    if (y > maxY) maxY = y;
	    if (y < minY) minY = y;

	    printf ("#=# %.12e %.12e\n", x, y);
	  }
	}
	
	//lb -= 1;
	//ub += 1;

	// plot lines defining constraint (only for cuts involving at
	// most two variables (that is, w is a unary function)
	for (int jj=j; jj < cs.sizeRowCuts (); jj++) {

	  CouNumber lb0 = lb, 
	            ub0 = ub;

	  const double *el  = cs.rowCutPtr (jj) -> row (). getElements ();
	  double  rhs = cs.rowCutPtr (jj) -> rhs ();

	  if (fabs (el [1]) > COUENNE_EPS) {
	    lb0 = mymax (lb, mymin ((rhs - el [0] * minY) / el [1], (rhs - el [0] * maxY) / el [1]));
	    ub0 = mymin (ub, mymax ((rhs - el [0] * minY) / el [1], (rhs - el [0] * maxY) / el [1]));
	  }

	  printf ("#=# #m=2,S=%d\n", (cs.rowCutPtr (jj) -> sense () == 'L') ? 10:11);

	  printf ("#=# %.12e %.12e\n", lb0, (rhs - el [1] * lb0) / el [0]);
	  printf ("#=# %.12e %.12e\n", ub0, (rhs - el [1] * ub0) / el [0]);
	}

	cg -> Problem () -> X () [xi] = curx;
      }
    }
  }
}

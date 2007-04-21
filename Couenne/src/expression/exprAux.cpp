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

#include <unistd.h> // for the sleep() in generateCuts()

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

  int j = cs.sizeRowCuts ();
  CouNumber l;

  if ((!(cg -> isFirst ())) && 
      (fabs ((l = expression::Ubound (varIndex_)) - 
	          expression::Lbound (varIndex_)) < COUENNE_EPS)) {

    cg -> createCut (cs, l, 0, varIndex_, 1.);
  } 
  else image_ -> generateCuts (this, si, cs, cg);

  //  if (!(cg -> isFirst ())) 
  if (j < cs.sizeRowCuts ())
  if (0)
    {
      printf ("----------------Generated cut for "); 
      print (std::cout);  printf (" := ");
      image_ -> print (std::cout); 
      printf("\n");
      for (int jj=j; jj < cs.sizeRowCuts ();jj++)
	cs.rowCutPtr (jj) -> print ();
    }


  if (0)  // [cool!] print graph-readable output for displaying
          // inequalities on a Cartesian plane

    if ((image_ -> code () == COU_EXPRSIN) || 
	(image_ -> code () == COU_EXPRCOS)) {

      printf (" ==> "); print (std::cout); printf ("\n");

      expression *lbe, *ube;

      int xi = image_ -> Argument () -> Index ();
      printf ("looking into w_%d = f (x_%d)\n", varIndex_, xi);

      image_ -> Argument () -> getBounds (lbe, ube);


      CouNumber lb   = (*lbe) (),
	        ub   = (*ube) ();

      delete lbe;
      delete ube;

      if (xi >= 0) {

	CouNumber curx = expression::Variable (xi);

#define N_STEPS 100

	// plot function
	for (CouNumber x = lb; 
	     x <= ub; 
	     x += ((ub - lb) / N_STEPS)) {

	  //	  expression::Variable (xi) = x;
	  printf ("#=# %.3f %.3f\n", x, (*image_) ());
	}

	// plot lines defining constraint (only for cuts involving at
	// most two variables (that is, w is a unary function)
	for (int jj=j; jj < cs.sizeRowCuts ();jj++) {
	  //	  const int    *ind = cs.rowCutPtr (jj) -> row (). getIndices  ();
	  const double *el  = cs.rowCutPtr (jj) -> row (). getElements ();
	  double  rhs = cs.rowCutPtr (jj) -> rhs ();

	  printf ("#=# #m=1,S=%d\n", (cs.rowCutPtr (jj) -> sense () == 'L') ? 10:11);

	  printf ("#=# %.3f %.3f\n", lb, (rhs - el [1] * lb) / el [0]);
	  printf ("#=# %.3f %.3f\n", ub, (rhs - el [1] * ub) / el [0]);
	}

	//	expression::Variable (xi) = curx;
      }
    }
}

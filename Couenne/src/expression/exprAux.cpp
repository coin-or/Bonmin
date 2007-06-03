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


/// I/O
void exprAux::print (std::ostream &out, bool descend, CouenneProblem *p) const {

  if (descend) 
    image_ -> print (out, descend, p);
  else out << "w_" << varIndex_;
}


// generate cuts for expression associated with this auxiliary

void exprAux::generateCuts (const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg, 
			    t_chg_bounds *chg) {

  int j = cs.sizeRowCuts ();
  CouNumber l, u;

  if ((!(cg -> isFirst ())) && 
      (fabs ((l = expression::Lbound (varIndex_)) - 
             (u = expression::Ubound (varIndex_))) < COUENNE_EPS))
    cg -> createCut (cs, (l+u)/2., 0, varIndex_, 1.);
  else image_ -> generateCuts (this, si, cs, cg, chg);

  // check if cuts have coefficients, rhs too large or too small

  for (int jj=j; jj < cs.sizeRowCuts (); jj++) {

    int           n   = cs.rowCutPtr (jj) -> row (). getNumElements();
    const double *el  = cs.rowCutPtr (jj) -> row (). getElements ();
    const int    *ind = cs.rowCutPtr (jj) -> row (). getIndices ();
    double        rhs = cs.rowCutPtr (jj) -> rhs ();

    while (n--)
      if      (fabs (el [n]) > COU_MAX_COEFF)  {
	printf ("coefficient too large: %.12e x%d\n", el [n], ind [n]);
      }
    //   else if (fabs (el[n]) < 1e-10) printf ("coefficient too small: %g x%d\n", el [n], ind [n]);

    if      (fabs (rhs) > COU_MAX_COEFF) {
      printf ("rhs too large: %.12e\n", rhs);
    }
    //    else if (fabs (rhs) < 1e-10) printf ("rhs too small: %g\n", rhs);
  }

  //  if (!(cg -> isFirst ())) 
  if (j < cs.sizeRowCuts ())
    if (0)
    {
      printf ("----------------Generated cut for "); 
      print (std::cout);  printf (" := ");
      image_ -> print (std::cout); 

      printf (" [%.3e,%.3e] <--- ", 
	      expression::Lbound (varIndex_), 
	      expression::Ubound (varIndex_));

      int index;
      if ((image_ -> Argument ()) && 
	  ((index = image_ -> Argument () -> Index ()) >= 0))
	printf ("[%.3e,%.3e] ", 
		expression::Lbound (index), 
		expression::Ubound (index));
      else if (image_ -> ArgList ())
	for (int i=0; i<image_ -> nArgs (); i++)
	  if ((index = image_ -> ArgList () [i] -> Index ()) >= 0)
	printf ("[%.3e,%.3e] ", 
		expression::Lbound (index), 
		expression::Ubound (index));
	
      printf("\n");
      for (int jj=j; jj < cs.sizeRowCuts ();jj++)
	cs.rowCutPtr (jj) -> print ();
    }

  //////////////////////////////////////////////////////////////

  if (0)
    draw_cuts (cs, cg, j, this, image_);
}

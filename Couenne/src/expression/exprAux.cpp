/*
 * Name:    exprAux.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the class of auxiliary variables
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneCutGenerator.hpp>
#include <CouenneTypes.h>
#include <expression.hpp>
#include <exprAux.hpp>
#include <exprVar.hpp>
#include <CouenneProblem.hpp>

//#define DEBUG


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
			    t_chg_bounds *chg, int,
			    CouNumber, CouNumber) {
#ifdef DEBUG
  static bool warned_large_coeff = false;

  int nrc = cs.sizeRowCuts (), ncc = cs.sizeColCuts ();
#endif

  /*
  if ((!(cg -> isFirst ())) && 
      ((l = expression::Lbound (varIndex_)) > -COUENNE_INFINITY) &&
      ((u = expression::Lbound (varIndex_)) <  COUENNE_INFINITY) &&
      (fabs (u-l) < COUENNE_EPS))
    cg -> createCut (cs, (l+u)/2., 0, varIndex_, 1.);
  else 
  */
  image_ -> generateCuts (this, si, cs, cg, chg);

  // check if cuts have coefficients, rhs too large or too small

#ifdef DEBUG
  if (warned_large_coeff)
    for (int jj=nrc; jj < cs.sizeRowCuts (); jj++) {

      OsiRowCut        *cut = cs.rowCutPtr (jj);
      CoinPackedVector  row = cut -> row ();

      int           n   = cut -> row (). getNumElements();
      const double *el  = row. getElements ();
      const int    *ind = row. getIndices ();
      double        rhs = cut -> rhs ();

      while (n--) {
	if (fabs (el [n]) > COU_MAX_COEFF)  {
	  printf ("### Warning: coefficient too large: %g x%d [", el [n], ind [n]);
	  cut -> print (); 
	  warned_large_coeff = true;
	  break;
	}

	if (fabs (rhs) > COU_MAX_COEFF) {
	  printf ("rhs too large %g: ", rhs);
	  cut -> print ();
	  warned_large_coeff = true;
	  break;
	}
      }
    }

  //  if (!(cg -> isFirst ())) 
  if ((nrc < cs.sizeRowCuts ()) || 
      (ncc < cs.sizeColCuts ()))
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
      for (int jj = nrc; jj < cs.sizeRowCuts (); jj++) cs.rowCutPtr (jj) -> print ();
      for (int jj = ncc; jj < cs.sizeColCuts (); jj++) cs.colCutPtr (jj) -> print ();
    }
#endif

  //////////////////////////////////////////////////////////////

#if 0
  draw_cuts (cs, cg, nrc, this, image_);
#endif
}

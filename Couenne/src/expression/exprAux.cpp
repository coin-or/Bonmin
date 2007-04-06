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

  int j = cs.sizeRowCuts ();
  CouNumber l; 

  if ((!(cg -> isFirst ())) && 
      (fabs ((l = expression::Ubound (varIndex_)) - 
	          expression::Lbound (varIndex_)) < COUENNE_EPS)) {

    OsiRowCut *cut = cg -> createCut (l, 0, varIndex_, 1.);
    if (cut) cs.insert (cut);
  } 
  else image_ -> generateCuts (this, si, cs, cg);

  //  if (!(cg -> isFirst ())) 
  //  if (j < cs.sizeRowCuts ())
  if (0)
    {
      printf ("----------------Generated cut for "); 
      print (std::cout);  printf (" := ");
      image_ -> print (std::cout); 
      printf("\n");
      for (;j < cs.sizeRowCuts ();j++)
	cs.rowCutPtr (j) -> print ();
    }
}

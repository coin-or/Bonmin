/*
 * Name:    problemIO.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the class CouenneProblem
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>
#include <fstream>
#include <iomanip> // to use the setprecision manipulator

#include <expression.hpp>
#include <exprAux.hpp>

#include <CouenneProblem.hpp>


// output content of the problem

void CouenneProblem::print (std::ostream &out) {

  printf ("objectives:\n");
  for (std::vector <CouenneObjective *>::iterator i = objectives_.begin ();
       i != objectives_.end (); i++)
    (*i) -> print (out);

  printf ("constraints:\n");
  for (std::vector <CouenneConstraint *>::iterator i = constraints_.begin ();
       i != constraints_.end (); i++)
    (*i) -> print (out);

  /*printf ("auxiliaries:\n");
  for (std::vector <exprAux *>::iterator i = auxiliaries_.begin ();
  i != auxiliaries_.end (); i++) */

  printf ("variables:\n");
  for (std::vector <exprVar *>::iterator i = variables_.begin ();
       i != variables_.end (); i++) 
    if ((*i) -> Type () == AUX) {
      if ((*i) -> Multiplicity () > 0) {

	(*i) -> print (out);
	out << " [" << (*i) -> rank (NULL) 
	    << ","  << (*i) -> Multiplicity () << "] := ";
	if ((*i) -> Image ())
	  (*i) -> Image () -> print (out, false, this); 
	out << " [ " << (*((*i) -> Lb ())) (); //-> print (out);
	out << " , " << (*((*i) -> Ub ())) (); //-> print (out);
	out << " ] " << std::endl;
      }
    } else {
      (*i) -> print (out);
      out << " in [" 
	  << lb_ [(*i) -> Index ()] << ',' 
	  << ub_ [(*i) -> Index ()] << "]\n";
    }

  if (optimum_) {
    printf ("best known solution: (%g", *optimum_);
    for (int i=1; i < nVars (); i++)
      printf (",%g", optimum_ [i]);
    printf (")\n");
  }

  if (fabs (bestObj_) < COUENNE_INFINITY)
    printf ("best known objective: %g\n", bestObj_);

  printf ("end\n");
} 


/// read optimal solution into member optimum
bool readOptimum (const std::string &fname, 
		  CouNumber *& optimum, 
		  CouNumber &bestObj, 
		  CouenneProblem *problem) {

  int nvars = problem -> nVars (),
      nOrig = problem -> nOrig ();
  //      nauxs = problem -> nAuxs ();

  FILE *f = fopen (fname.c_str (), "r");
  if (!f) return false;

  optimum = (CouNumber *) realloc (optimum, nvars * sizeof (CouNumber));

  if (fscanf (f, "%lf", &bestObj) < 1) 
    return false;

  for (int i = 0; i < nOrig; i++)
    if (fscanf (f, "%lf", optimum + i) < 1) 
      return false;

  // save current expression vectors
  CouNumber
    *xS = expression::Variables (),
    *lS = expression::Lbounds   (),
    *uS = expression::Ubounds   ();

  // now propagate value to aux variables
  expression::update (optimum, NULL, NULL);

  // only one loop is sufficient here, since auxiliary variables are
  // defined in such a way that w_i does NOT depend on w_j if i<j.

  for (register int i = 0, j = problem -> nVars (); j--; i++)
    if (problem -> Var (i) -> Type () == AUX)
      optimum [problem -> Var (i) -> Index ()] = (*(problem -> Var (i) -> Image ())) ();

  expression::update (problem -> X  (), 
		      problem -> Lb (), 
		      problem -> Ub ());

  // restore previous value/bound vectors
  expression::update (xS, lS, uS);

  return true;
}

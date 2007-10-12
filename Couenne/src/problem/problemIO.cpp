/*
 * Name:    problemIO.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the class CouenneProblem
 *
 * (C) Carnegie-Mellon University, 2006. 
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
       i != objectives_.end (); ++i)
    (*i) -> print (out);

  printf ("constraints:\n");
  for (std::vector <CouenneConstraint *>::iterator i = constraints_.begin ();
       i != constraints_.end (); ++i)
    (*i) -> print (out);

  printf ("variables:\n");
  for (std::vector <exprVar *>::iterator i = variables_.begin ();
       i != variables_.end (); ++i) 
    if ((*i) -> Type () == AUX) {
      if ((*i) -> Multiplicity () > 0) {

	(*i) -> print (out);
	out << " [" << (*i) -> rank (NULL) 
	    << ","  << (*i) -> Multiplicity () << "] := ";
	if ((*i) -> Image ())
	  (*i) -> Image () -> print (out, false, this); 
	out << " [ " << (*((*i) -> Lb ())) ();
	out << " , " << (*((*i) -> Ub ())) ();
	out << " ] " << std::endl;
	/*
	expression *lb, *ub;
	(*i) -> getBounds (lb, ub);
	out << " {";  lb -> print (out);
	out << " , "; ub -> print (out);
	out << "}\n"; 
	*/
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

  // TODO: this procedure is crippled by the new auxiliary handling
  // which replaces original variables with auxiliaries. The problem
  // could be that some originally auxiliary variables (er) do not
  // have an optimal value but are evaluated as independent.

  // Actually, forget the above. It can only happen in extended
  // formulations, whose original variables are either original or
  // auxiliary in the original problem.

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

  for (register int i = 0, j = problem -> nVars (); j--; i++) {
    exprVar *var = problem -> Var (problem -> evalOrder (i));
    if (var -> Type () == AUX)
      optimum [var -> Index ()] = (*(var -> Image ())) ();
  }

  // restore previous value/bound vectors
  expression::update (xS, lS, uS);

  return true;
}

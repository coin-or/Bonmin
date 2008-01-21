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

#include "expression.hpp"
#include "exprAux.hpp"

#include "CouenneProblem.hpp"


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

    if (((*i) -> Type () != AUX) || ((*i) -> Multiplicity () > 0)) {

      (*i) -> print (out);

      if (((*i) -> Type () == AUX) && ((*i) -> Multiplicity () > 0)) {

	//printf (" %x", (*i));
	out << " (r:" << (*i) -> rank () 
	    << ", m:"  << (*i) -> Multiplicity () << ") := ";
	if ((*i) -> Image ())
	  (*i) -> Image () -> print (out, false); 
      }

      if ((fabs ((*((*i) -> Lb ())) ())     < COUENNE_EPS) &&
	  (fabs ((*((*i) -> Ub ())) () - 1) < COUENNE_EPS) &&
	  (*i) -> isInteger ()) out << " binary";
      else {

	out << " [ " << (*((*i) -> Lb ())) ();
	out << " , " << (*((*i) -> Ub ())) ();
	out << " ]";
	/*
	if ((*i) -> Image ()) {
	  expression *lb, *ub;
	  (*i) -> Image () -> getBounds (lb, ub);
	  out << " {";  lb -> print (out);
	  out << " , "; ub -> print (out);
	  out << "} "; 
	}
	*/
	if ((*i) -> isInteger ()) out << " integer";
      }

      out << std::endl;
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
bool CouenneProblem::readOptimum (const std::string &fname) {

  // TODO: this procedure is crippled by the new auxiliary handling
  // which replaces original variables with auxiliaries. The problem
  // could be that some originally auxiliary variables (er...) do not
  // have an optimal value but are evaluated as independent.

  // Actually, forget the above. It can only happen in extended
  // formulations, whose original variables are either original or
  // auxiliary in the original problem.

  int nvars = nVars ();

  FILE *f = fopen (fname.c_str (), "r");
  if (!f) return false;

  optimum_ = (CouNumber *) realloc (optimum_, nvars * sizeof (CouNumber));

  CoinFillN (optimum_, nvars, 0.);

  if (fscanf (f, "%lf", &bestObj_) < 1) 
    return false;

  for (int i = 0; i < nOrig_; i++)
    if (fscanf (f, "%lf", optimum_ + i) < 1) 
      return false;

  // save current expression vectors
  CouNumber
    *xS = expression::Variables (),
    *lS = expression::Lbounds   (),
    *uS = expression::Ubounds   ();

  // now propagate value to aux variables
  expression::update (optimum_, NULL, NULL);

  // only one loop is sufficient here, since auxiliary variables are
  // defined in such a way that w_i does NOT depend on w_j if i<j.

  for (int i = 0, j = nVars (); j--; i++) {
    exprVar *var = variables_ [numbering_ [i]];
    if (var -> Type () == AUX)
      optimum_ [var -> Index ()] = (*(var -> Image ())) ();
  }

  // restore previous value/bound vectors
  expression::update (xS, lS, uS);

  return true;
}


/// read cutoff into member optimum
void CouenneProblem::readCutoff (const std::string &fname) {

  CouNumber val;

  FILE *f = fopen (fname.c_str (), "r");
  if (!f) return;

  if (fscanf (f, "%lf", &val) < 1) 
    return;

  setCutOff (val);
}

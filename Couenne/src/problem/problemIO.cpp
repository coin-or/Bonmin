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

#include <CouenneTypes.h>
#include <expression.h>
#include <exprConst.h>
#include <exprGroup.h>
#include <exprAux.h>

#include <CouenneProblem.h>
#include <CouenneProblemElem.h>


// output content of the problem

void CouenneProblem::print (std::ostream &out = std::cout) {

  printf ("objectives:\n");
  for (std::vector <Objective *>::iterator i = objectives_.begin ();
       i != objectives_.end (); i++)
    (*i) -> print (out);

  printf ("constraints:\n");
  for (std::vector <CouenneConstraint *>::iterator i = constraints_.begin ();
       i != constraints_.end (); i++)
    (*i) -> print (out);

  printf ("auxiliaries:\n");
  for (std::vector <exprAux *>::iterator i = auxiliaries_.begin ();
       i != auxiliaries_.end (); i++) 

    if ((*i) -> Multiplicity () > 0) {

      (*i) -> print (out);
      out << " [" << (*i) -> rank (NULL) 
	  << ","  << (*i) -> Multiplicity () << "] := ";
      (*i) -> Image () -> print (out, false, this); 
      out << " [ " << (*((*i) -> Lb ())) (); //-> print (out);
      out << " , " << (*((*i) -> Ub ())) (); //-> print (out);
      out << " ] " << std::endl;
    }

  printf ("bounds:\n");
  for (std::vector <exprVar *>::iterator i = variables_.begin ();
       i != variables_.end (); i++) {
    (*i) -> print (out);
    out << " in [" 
	<< lb_ [(*i) -> Index ()] << ',' 
	<< ub_ [(*i) -> Index ()] << "]\n";
  }

  if (optimum_) {
    printf ("best known solution: (%g", *optimum_);
    for (int i=1; i < nVars () + nAuxs (); i++)
      printf (",%g", optimum_ [i]);
    printf (")\n");
  }

  if (fabs (bestObj_) < COUENNE_INFINITY)
    printf ("best known objective: %g\n", bestObj_);

  printf ("end\n");
} 


// store problem in a .mod file (AMPL)

void CouenneProblem::writeMod (const std::string &fname,  /// name of the mod file
			       bool aux) {                /// with or without auxiliaries?

  std::ofstream f (fname.c_str ());

  f << std::setprecision (10);

  // original variables, integer and non //////////////////////////////////////////////

  f << "# Problem name: " << fname << std::endl << std::endl 
    << "# original variables" << std::endl << std::endl;

  for (int i=0; i < nVars (); i++) {

    f << "var ";
    variables_ [i] -> print (f);
    if (lb_ [i] > - COUENNE_INFINITY + 1) f << " >= " << lb_ [i];
    if (ub_ [i] < + COUENNE_INFINITY - 1) f << " <= " << ub_ [i];
    if (variables_ [i] -> isInteger ())   f << " integer";
    f << " default " << x_ [i] << ';' << std::endl;
  }


  // defined (aux) variables, declaration ///////////////////////////////////////////

  if (aux) {

    initAuxs (x_, lb_, ub_);

    f << std::endl << "# auxiliary variables" << std::endl << std::endl;

    for (std::vector <exprAux *>::iterator i = auxiliaries_.begin ();
	 i != auxiliaries_.end ();
	 i++) 

      if ((*i) -> Multiplicity () > 0) {

	f << "var "; (*i) -> print (f, false, this);
	//    f << " = ";  (*i) -> Image () -> print (f);
	CouNumber bound;

	if ((bound = (*((*i) -> Lb ())) ()) > - COUENNE_INFINITY) f << " >= " << bound;
	if ((bound = (*((*i) -> Ub ())) ()) <   COUENNE_INFINITY) f << " <= " << bound;
	if ((*i) -> isInteger ()) f << " integer";
	f << " default " << (*((*i) -> Image ())) () << ';' << std::endl;
      }
  }


  // objective function /////////////////////////////////////////////////////////////

  f << std::endl << "# objective" << std::endl << std::endl;

  if (objectives_ [0] -> Sense () == MINIMIZE) f << "minimize";
  else                                         f << "maximize";

  f << " obj: ";  
  objectives_ [0] -> Body () -> print (f, !aux, this); 
  f << ';' << std::endl; 


  // defined (aux) variables, with formula ///////////////////////////////////////////

  if (aux) {
    f << std::endl << "# aux. variables defined" << std::endl << std::endl;

    for (int i=0; i < nAuxs (); i++) 
      if (auxiliaries_ [i] -> Multiplicity () > 0) {

	f << "aux" << i << ": "; auxiliaries_ [i] -> print (f, false, this);
	f << " = ";  

	auxiliaries_ [i] -> Image () -> print (f, false, this);
	f << ';' << std::endl;
      }
  }


  // write constraints //////////////////////////////////////////////////////////////

  f << std::endl << "# constraints" << std::endl << std::endl;

  if (!aux) // print h_i(x,y) <= ub, >= lb
    for (std::vector <exprAux *>::iterator i = auxiliaries_.begin ();
	 i != auxiliaries_.end ();
	 i++) 

      if ((*i) -> Multiplicity () > 0) {
	
	CouNumber bound;

	if ((bound = (*((*i) -> Lb ())) ()) > - COUENNE_INFINITY) {
	  f << "conAuxLb" << (*i) -> Index () << ": ";
	  (*i) -> print (f, true, this);
	  f << ">= " << bound << ';' << std::endl;
	}

	if ((bound = (*((*i) -> Ub ())) ()) <   COUENNE_INFINITY) {
	  f << "conAuxUb" << (*i) -> Index () << ": ";
	  (*i) -> print (f, true, this);
	  f << "<= " << bound << ';' << std::endl;
	}
      }


  for (int i=0; i < nNLCons (); i++) {

    // get numerical value of lower, upper bound
    CouNumber lb = (constraints_ [i] -> Lb () -> Value ()),
              ub = (constraints_ [i] -> Ub () -> Value ());

    f << "con" << i << ": ";
    constraints_ [i] -> Body () -> print (f, !aux, this);

    if (lb > - COUENNE_INFINITY + 1) {
      f << ' ';
      if (fabs (ub-lb) > COUENNE_EPS) 
	f << '>';
      f << "= " << lb << ';' << std::endl;
    }
    else  f << " <= " << ub << ';' << std::endl;

    // if range constraint, print it once again

    if ((   lb > - COUENNE_INFINITY + 1) 
	&& (ub <   COUENNE_INFINITY - 1)
	&& (fabs (ub-lb) > COUENNE_EPS)) {

      f << "con" << i << "_rng: ";
      constraints_ [i] -> Body () -> print (f, !aux, this);
      f << " <= " << ub << ';' << std::endl;
    }
  }
}

/// read optimal solution into member optimum
bool readOptimum (const std::string &fname, 
		  CouNumber *& optimum, 
		  CouNumber &bestObj, 
		  CouenneProblem *problem) {

  int nvars = problem -> nVars (),
      nauxs = problem -> nAuxs ();

  FILE *f = fopen (fname.c_str (), "r");
  if (!f) return false;

  optimum = (CouNumber *) realloc (optimum, (nvars + nauxs) * sizeof (CouNumber));

  if (fscanf (f, "%lf", &bestObj) < 1) 
    return false;

  for (int i=0; i<nvars; i++)
    if (fscanf (f, "%lf", optimum + i) < 1) 
      return false;

  // now propagate value to aux variables
  expression::update (optimum, NULL, NULL);

  // only one loop is sufficient here, since auxiliary variables are
  // defined in such a way that w_i does NOT depend on w_j if i<j.

  for (register int i = 0, j = problem -> nVars (); i < problem -> nAuxs (); i++, j++)
    optimum [j] = (*(problem -> Aux (i) -> Image ())) ();

  expression::update (problem -> X  (), 
		      problem -> Lb (), 
		      problem -> Ub ());

  return true;
}

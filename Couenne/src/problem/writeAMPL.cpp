/*
 * Name:    extended.cpp
 * Author:  Pietro Belotti
 * Purpose: save problem in AMPL format
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <fstream>
#include <iomanip> // to use the setprecision manipulator

#include "CouenneProblem.hpp"


// store problem in a .mod file (AMPL)

void CouenneProblem::writeAMPL (const std::string &fname,  /// name of the mod file
				bool aux) {                /// with or without auxiliaries?

  std::ofstream f (fname.c_str ());

  f << std::setprecision (10);

  // original variables, integer and non //////////////////////////////////////////////

  f << "# Problem name: " << fname << std::endl << std::endl 
    << "# original variables" << std::endl << std::endl;

  for (int i=0; i < nVars (); i++) {

    f << "var ";
    variables_ [i] -> print (f);
    if (Lb (i) > - COUENNE_INFINITY + 1) f << " >= " << Lb (i);
    if (Ub (i) < + COUENNE_INFINITY - 1) f << " <= " << Ub (i);
    if (variables_ [i] -> isInteger ())   f << " integer";
    if (fabs (X (i)) < COUENNE_INFINITY)    
      f << " default " << X (i); 
    f << ';' << std::endl;
  }


  // defined (aux) variables, declaration ///////////////////////////////////////////
  /*
  if (aux) {

    initAuxs (x_, lb_, ub_);

    f << std::endl << "# auxiliary variables" << std::endl << std::endl;

    for (std::vector <exprVar *>::iterator i = variables_.begin ();
	 i != variables_.end ();
	 i++) 

      if ((*i) -> Type () == AUX) {

	exprAux *aux = dynamic_cast <exprAux *> (*i);

	if (aux -> Multiplicity () > 0) {

	  f << "var "; (*i) -> print (f, false, this);
	  //    f << " = ";  (*i) -> Image () -> print (f);
	  CouNumber bound;

	  if ((bound = (*((*i) -> Lb ())) ()) > - COUENNE_INFINITY) f << " >= " << bound;
	  if ((bound = (*((*i) -> Ub ())) ()) <   COUENNE_INFINITY) f << " <= " << bound;
	  if ((*i) -> isInteger ()) f << " integer";

	  f << " default " << (*((*i) -> Image ())) () << ';' << std::endl;
	}
      }
  }
  */

  // objective function /////////////////////////////////////////////////////////////

  f << std::endl << "# objective" << std::endl << std::endl;

  if (objectives_ [0] -> Sense () == MINIMIZE) f << "minimize";
  else                                         f << "maximize";

  f << " obj: ";  
  objectives_ [0] -> Body () -> print (f, !aux); 
  f << ';' << std::endl; 


  // defined (aux) variables, with formula ///////////////////////////////////////////

  if (aux) {

    f << std::endl << "# aux. variables defined" << std::endl << std::endl;

    for (int i=0; i < nVars (); i++)

      if ((variables_ [i] -> Type () == AUX) && 
	  (variables_ [i] -> Multiplicity () > 0)) {

	f << "aux" << i << ": "; variables_ [i] -> print (f, false);
	f << " = ";  

	variables_ [i] -> Image () -> print (f, false);
	f << ';' << std::endl;
      }
  }


  // write constraints //////////////////////////////////////////////////////////////

  f << std::endl << "# constraints" << std::endl << std::endl;

  if (!aux) // print h_i(x,y) <= ub, >= lb
    for (std::vector <exprVar *>::iterator i = variables_.begin ();
	 i != variables_.end ();
	 ++i) 

      if (((*i) -> Type () == AUX) && 
	  ((*i) -> Multiplicity () > 0)) {
	
	CouNumber bound;

	if ((bound = (*((*i) -> Lb ())) ()) > - COUENNE_INFINITY) {
	  f << "conAuxLb" << (*i) -> Index () << ": ";
	  (*i) -> print (f, true);
	  f << ">= " << bound << ';' << std::endl;
	}

	if ((bound = (*((*i) -> Ub ())) ()) <   COUENNE_INFINITY) {
	  f << "conAuxUb" << (*i) -> Index () << ": ";
	  (*i) -> print (f, true);
	  f << "<= " << bound << ';' << std::endl;
	}
      }


  for (int i=0; i < nCons (); i++) {

    // get numerical value of lower, upper bound
    CouNumber lb = (constraints_ [i] -> Lb () -> Value ()),
              ub = (constraints_ [i] -> Ub () -> Value ());

    f << "con" << i << ": ";
    constraints_ [i] -> Body () -> print (f, !aux);

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
      constraints_ [i] -> Body () -> print (f, !aux);
      f << " <= " << ub << ';' << std::endl;
    }
  }
}

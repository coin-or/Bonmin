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
       i != auxiliaries_.end (); i++) {
    (*i) -> print (out);
    out << " [" << (*i) -> rank (NULL) 
	<< ","  << (*i) -> Multiplicity () << "] := ";
    (*i) -> Image () -> print (out, false, this); 
    out << " [ " << (*((*i) -> Lb ())) (); //-> print (out);
    out << " , " << (*((*i) -> Ub ())) (); //-> print (out);
    out << " ] " << std::endl;
  }
  printf ("end\n");
} 


// store problem in a .mod file (AMPL)

void CouenneProblem::writeMod (char *fname,  /// name of the mod file
			       bool aux) {   /// with or without auxiliaries?

  std::ofstream f (fname);

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
    f << ';' << std::endl;
  }


  // defined (aux) variables, declaration ///////////////////////////////////////////

  if (aux) {

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
	f << ';' << std::endl;
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

	//	if (auxiliaries_ [i] -> Image () -> code () == COU_EXPRGROUP)
	  //	  dynamic_cast <exprGroup *> (auxiliaries_ [i] -> Image ()) -> print (f, this);
	  //	else 
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
	
	//	f << "conAux"; (*i) -> print (f);
	//    f << " = ";  (*i) -> Image () -> print (f);
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

	//	if ((bound = (*((*i) -> Ub ())) ()) <   COUENNE_INFINITY) f << " <= " << bound;
	//	if ((*i) -> isInteger ()) f << " integer";
	//	f << ';' << std::endl;
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

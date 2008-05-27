/*
 * Name:    writeGAMS.cpp
 * Author:  Pietro Belotti
 * Purpose: save a problem in GAMS format
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <fstream>
#include <iomanip> // to use the setprecision manipulator

#include "CouenneProblem.hpp"


/// Write nonlinear problem to a .gms file
/// 
/// @param fname Name of the .gms file to be written.
void CouenneProblem::writeGAMS (const std::string &fname) {

  std::ofstream f (fname.c_str ());

  f << std::setprecision (10);

  const int nline = 20;

  // header //////////////////////////////////////////////////////////////

  f << "* MINLP Written by Couenne (https://projects.coin-or.org/Bonmin/wiki/BonCouenne)" << std::endl
    << "* " << nVars () << " variables, " << nCons () << " constraints" << std::endl << std::endl;

  // preamble ////////////////////////////////////////////////////////////

  f << "  variables   ";

  for (int i=0, j=0; i < nVars (); i++) 
    if (Var (i) -> Multiplicity () > 0) {

      if (j) f << ','; 
      if (!(++j % nline)) f << std::endl << "             ";

      Var (i) -> print (f);

      /*if (Var (i) -> isInteger ())
	if ((fabs (Var (i) -> lb ())     < COUENNE_EPS) && 
	    (fabs (Var (i) -> ub () - 1) < COUENNE_EPS))
	  f << 'b';
	else f << 'i';
      else f << 'x';
      f << Var (i) -> Index ();*/
    }

  f << ",objvar;" << std::endl << std::endl;

  bool 
    have_positive = false, 
    have_binary   = false,
    have_integer  = false;

  // check if there are positive variables
  for (int i=0; i < nVars (); i++)
    if ((Var (i) -> Multiplicity () > 0) &&
	(fabs (Var (i) -> lb ()) < COUENNE_EPS)  &&
	!(Var (i) -> isInteger () &&
	  (fabs (Var (i) -> ub () - 1) < COUENNE_EPS))) {

      have_positive = true;
      break;
    }

  // check if there are binary variables
  for (int i=0; i < nVars (); i++)
    if ((Var (i) -> Multiplicity () > 0) &&
	Var (i) -> isInteger () &&
	(fabs (Var (i) -> lb ())     < COUENNE_EPS) && 
	(fabs (Var (i) -> ub () - 1) < COUENNE_EPS)) {

      have_binary = true;
      break;
    }

  // check if there are integer variables
  for (int i=0; i < nVars (); i++)
    if ((Var (i) -> Multiplicity () > 0) &&
	(Var (i) -> isInteger ()) &&
	((fabs (Var (i) -> lb ())     > COUENNE_EPS) ||
	 (fabs (Var (i) -> ub () - 1) > COUENNE_EPS))) {

      have_integer = true;
      break;
    }

  // positive /////////////////////////////////////////////////////////////////
  if (have_positive) {

    f << "  positive variables "; 
    // those with lower bound exactly zero (to save space...)

    for (int i=0, j=0; i < nVars (); i++) 
      if ((Var (i) -> Multiplicity () > 0) && 
	  (fabs (Var (i) -> lb ()) < COUENNE_EPS)  &&
	  !(Var (i) -> isInteger () &&
	    (fabs (Var (i) -> ub () - 1) < COUENNE_EPS))) {

	if (j) f << ','; 
	if (!(++j % nline)) f << std::endl << "             ";

	Var (i) -> print (f);

	/*if (Var (i) -> isInteger ()) f << 'i';
	else                         f << 'x';
	f << Var (i) -> Index ();*/
      }

    f << ';' << std::endl << std::endl;
  }

  // binary /////////////////////////////////////////////////////////////////
  if (have_binary) {

    f << "  binary variables ";

    for (int i=0, j=0; i < nVars (); i++) 
      if ((Var (i) -> Multiplicity () > 0) && 
	  Var (i) -> isInteger () &&
	  (fabs (Var (i) -> lb ())     < COUENNE_EPS) && 
	  (fabs (Var (i) -> ub () - 1) < COUENNE_EPS)) {

	if (j) f << ','; 
	if (!(++j % nline)) f << std::endl << "             ";
	Var (i) -> print (f);
	//f << 'b' << Var (i) -> Index ();
      }
    f << ';' << std::endl << std::endl;
  }

  // integer /////////////////////////////////////////////////////////////////
  if (have_integer) {

    f << "  integer variables ";

    for (int i=0, j=0; i < nVars (); i++) 
      if ((Var (i) -> Multiplicity () > 0) &&
	  (Var (i) -> isInteger ()) &&
	  ((fabs (Var (i) -> lb ())     > COUENNE_EPS) ||
	   (fabs (Var (i) -> ub () - 1) > COUENNE_EPS))) {

	if (j) f << ','; 
	if (!(++j % nline)) f << std::endl << "             ";
	Var (i) -> print (f);
	//f << 'i' << Var (i) -> Index ();
      }
    f << ';' << std::endl << std::endl;
  }

  f << "  equations ";

  int i=0, j=0;

  bool no_eqns = true;

  for (; i < nCons (); i++) {

    if (j) f << ','; 
    if (!(++j % nline)) f << std::endl << "                  ";
    f << 'e' << j;
    no_eqns = false;
  }

  for (; i < nVars (); i++) 
    if ((Var (i) -> Type () == AUX) &&
	(Var (i) -> Multiplicity () > 0)) {

      if (j) f << ','; 
      if (!(++j % nline)) f << std::endl << "                ";
      f << 'e' << j;
      no_eqns = false;
    }

  if (!no_eqns) f << ',';
  f << "objdef;" << std::endl << std::endl;

  // equations ///////////////////////////////////////////////////////////

  i=j=0;

  for (; i < nCons (); i++) {

    f << 'e' << ++j << "..  ";

    CouNumber 
      lb = (*(Con (i) -> Lb ())) (),
      ub = (*(Con (i) -> Ub ())) ();

    assert ((lb > -COUENNE_INFINITY) || 
	    (ub <  COUENNE_INFINITY));

    if (fabs (lb - ub) < COUENNE_EPS) {
      Con (i) -> Body () -> print (f);
      f << " =E= " << lb;
    } else if (lb > -COUENNE_INFINITY) {
      Con (i) -> Body () -> print (f);
      f << " =G= " << lb;
    } else if (ub <  COUENNE_INFINITY) {
      Con (i) -> Body () -> print (f);
      f << " =L= " << ub;
    }

    f << ';' << std::endl << std::endl;
  }

  for (; i < nVars (); i++) 
    if ((Var (i) -> Type () == AUX) &&
	(Var (i) -> Multiplicity () > 0)) {

      f << 'e' << ++j << "..  ";
      f << " - ";
      Var (i) -> print (f); 
      f << " + ";
      Var (i) -> Image () -> print ();
      f << " =E= 0;" << std::endl << std::endl;
    }

  f << "objdef.. - objvar + ";
  Obj (0) -> Body () -> print (f);
  f << " =E= 0;" << std::endl << std::endl;

  // initial levels //////////////////////////////////////////////////////

  for (int i=0; i < nVars (); i++) 
    if (Var (i) -> Multiplicity () > 0) {

      Var (i) -> print (f); 
      f << ".l = " << X (i) << ';' << std::endl;

      if ((fabs (Lb (i)) > COUENNE_EPS) && (Lb (i) > -COUENNE_INFINITY)) {
	Var (i) -> print (f); 
	f << ".lo = " << Lb (i) << ';' << std::endl;
      }

      if (Ub (i) < COUENNE_INFINITY) {
	Var (i) -> print (f); 
	f << ".up = " << Ub (i) << ';' << std::endl;
      }
    }

  // bounds //////////////////////////////////////////////////////////////

  f << "Model m / all /;" << std::endl
    << "m.limrow=0; m.limcol=0;" << std::endl
    << "Solve m using MINLP minimizing objvar;" << std::endl;

  f.close ();
}

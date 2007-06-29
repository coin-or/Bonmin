/*
 * Name:    updateSol.cpp
 * Author:  Pietro Belotti
 * Purpose: update best quadratic solution if necessary
 *
 *
 * (C) Pietro Belotti, Carnegie Mellon University. This file is
 * distributed under the Common Public License.
 */

#include <stdio.h>
#include <stdlib.h>
#include <OsiCpxSolverInterface.hpp>
#include <SdpCutGen.hpp>


void SdpCutGen::updateSol (OsiSolverInterface &si) {

  static int ncalls = 0;

  /*
    insert best solution of 30-60-1
  */

  double *sol = const_cast <double *> (si.getColSolution ());

#if 0
  if (ncalls < 20) {
    sol = new double [30];

    double s [30] = {0,0,1,0,0,0,1,1,1,0,0,1,1,1,1,0,1,0,0,1,0,1,0,1,1,0,1,1,1,0};

    for (int i=0; i<30; i++) 
      sol [i] = s [i];
  }
#endif

  double obj = 0;

  for (int i=0; i<n_; i++)
    obj += b_ [i] * sol [i];

  obj *= 2;

  for (int i=0; i<n_; i++)
    for (int j=0; j < n_; j++)
      obj += Q_ [i] [j] * sol [i] * sol [j];

  obj /= 2;

  currObj_ = obj;

  if (obj > bestObj_) {

    bestObj_ = obj;

    if (!bestSol_)
      bestSol_ = new double [n_];

    printf ("found new %.6f: ", bestObj_);

    for (int i=0; i<n_; i++)
      printf ("%.1g ", sol [i]);

    printf ("\n");

    for (int i=0; i<n_; i++)
      bestSol_ [i] = sol [i];
  }

  ncalls++;
}

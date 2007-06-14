/*
 * Apply semidefinite cuts to solve Unconstrained Quadratic
 * Programming problems
 *
 * (C) Pietro Belotti, Carnegie Mellon University
 */

#include <stdio.h>
#include <OsiCuts.hpp>
#include <OsiSolverInterface.hpp>

#define indexQ(i,j,n) ((n) + (i) * (2*(n)-1-(i)) / 2 + (j))

int createCut (OsiCuts &, double, int, 
	       int=-1, double=0,
	       int=-1, double=0,
	       int=-1, double=0, bool=false);

// fill in initial RLT
void populateProblem (int n, double *b, double **Q, OsiSolverInterface *si) {

  int N = n*(n+3)/2;

  // add all n + n(n+1)/2 variables
  for (register int i=N; i--;)
    si -> addCol (0, NULL, NULL, 0, 1, 0);

  OsiCuts cs;

  // initial square cuts:
  for (int i=0; i<n; i++) {
    int ind = indexQ (i,i,n);
    createCut (cs,  0, -1, ind, 1, i, -1); // upper segment  Xii <= xi 
    createCut (cs, -1, +1, ind, 1, i, -2); // lower envelope Xii >= 2xi - 1
  }

  // McCormick cuts
  for (int i=0; i<n; i++)
    for (int j=i+1; j<n; j++) {
      
      int ind = indexQ (i,j,n);
      createCut (cs, -1, +1, ind, 1, i, -1, j, -1); // Xij >= uj xi + ui xj - ui uj
      createCut (cs,  0, -1, ind, 1, i, -1);        // Xij <= uj xi + li xj - li uj
      createCut (cs,  0, -1, ind, 1, j, -1);        // Xij <= lj xi + ui xj - ui lj
    }

  // add constraints
  si -> applyCuts (cs);

  // create objective function
  double * obj = new double [N];

  for (int i=0; i<n; i++)
    obj [i] = b [i];

  for (int i=0; i<n; i++) {
    obj [indexQ (i,i,n)] = 0.5 * Q [i] [i];
    for (int j=i+1; j<n; j++)
      obj [indexQ (i,j,n)] = 0.5 * (Q [i] [j] + Q [j] [i]);
  }

  si -> setObjective (obj);
  si -> setObjSense  (-1);  // maximization

  delete [] obj;
}

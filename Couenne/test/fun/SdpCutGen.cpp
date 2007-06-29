#include <SdpCutGen.hpp>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <sys/time.h>
#include <sys/resource.h>

// sdpcut separator
void SdpCutGen::generateCuts (const OsiSolverInterface &si, OsiCuts &cs, 
			      const CglTreeInfo info) const {
  int ncuts;

  violated_ = false;

  printf ("================= Separation:\n");

#if 1
  static bool firstcall = true;
  if (firstcall) {
    ncuts = cs.sizeRowCuts ();
    separateGrad (si, cs);
    printf ("================= Grad: %3d cuts\n", cs.sizeRowCuts () - ncuts);
    firstcall = false;
  }
#endif

  ncuts = cs.sizeRowCuts ();
  separateTRM (si, cs);
  printf ("================= TRM: %3d cuts\n", cs.sizeRowCuts () - ncuts);

  ncuts = cs.sizeRowCuts ();
  separateEV  (si, cs);
  printf ("================= EV:  %3d cuts\n", cs.sizeRowCuts () - ncuts);

  //  ncuts = cs.sizeRowCuts ();
  //  separateLU (si, cs);
  //  printf ("================= LU:  %3d cuts\n", cs.sizeRowCuts () - ncuts);

  //  ncuts = cs.sizeRowCuts ();
  //  separateBK (si, cs);
  //  printf ("================= BK:  %3d cuts\n", cs.sizeRowCuts () - ncuts);
}

/// constructor
SdpCutGen::SdpCutGen  (int n, double *b, double **Q):
  n_ (n),
  currObj_ (-DBL_MAX),
  bestObj_ (-DBL_MAX),
  bestSol_ (NULL) {

  if (n <= 0) n = 1;

  b_ = new double   [n];
  Q_ = new double * [n];
  for (int i=n; i--;)
    Q_ [i] = new double [n];

  for (int i=n; i--;) {
    b_ [i] = b [i];
    for (int j=n; j--;) 
      Q_ [i] [j] = Q [i] [j];
  }
}

/// copy constructor
SdpCutGen::SdpCutGen  (const SdpCutGen &rhs):
  n_ (rhs.n_) {

  b_ = new double   [n_];
  Q_ = new double * [n_];
  for (int i=n_; i--;)
    Q_ [i] = new double [n_];

  for (int i=n_; i--;) {
    b_ [i] = rhs.b_ [i];
    for (int j=n_; j--;) 
      Q_ [i] [j] = rhs.Q_ [i] [j];
  }
}

/// destructor
SdpCutGen::~SdpCutGen () {

  while (n_--) delete [] Q_ [n_];
  delete [] Q_;
  delete [] b_; 
  if (bestSol_) delete [] bestSol_;
}

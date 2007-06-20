#include <SdpCutGen.hpp>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <sys/time.h>
#include <sys/resource.h>

// sdpcut separator
void SdpCutGen::generateCuts (const OsiSolverInterface &si, 
			      OsiCuts &cs, 
			      const CglTreeInfo info) const {
  static int ncalls = 0;

  int ncuts = cs.sizeRowCuts ();
  printf ("================= Separation:\n");
  separateTRM (si, cs);
  printf ("================= TRM: %3d cuts\n", cs.sizeRowCuts () - ncuts);
  ncuts = cs.sizeRowCuts ();
  separateEV  (si, cs);
  printf ("================= EV:  %3d cuts\n", cs.sizeRowCuts () - ncuts);
  ncuts = cs.sizeRowCuts ();
  //  separateLU (si, cs);
  printf ("================= LU:  %3d cuts\n", cs.sizeRowCuts () - ncuts);
  ncuts = cs.sizeRowCuts ();
  //  separateBK (si, cs);
  printf ("================= BK:  %3d cuts\n", cs.sizeRowCuts () - ncuts);

  ncalls++;
}

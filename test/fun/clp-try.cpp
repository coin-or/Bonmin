#include <OsiClpSolverInterface.hpp>

int main (int argc, char **argv) {

  OsiClpSolverInterface si;

  si.readLp (argv [1]);

  int nvars = si.getNumCols ();

  double *objcoe = new double [nvars];

  for (int i=0; i<nvars; i++)
    objcoe [i] = 0;

  si.initialSolve();

  for (int i=0; i<nvars; i++) {
    si.setObjSense (1);
    objcoe [i] = 1;
    si.setObjective (objcoe);
    si.resolve();
    if (si.isProvenOptimal())
      printf ("------------------------------min x%d: %g\n", i, si. getObjValue ());
    else printf (".............\n");
    objcoe [i] = 0;
  }

  printf ("#################################\n");

  for (int i=0; i<nvars; i++) {
    si.setObjSense (1);
    objcoe [i] = -1;
    si.setObjective (objcoe);
    si.resolve();
    if (si.isProvenOptimal())
      printf ("------------------------------max x%d: %g\n", i, - si. getObjValue ());
    else printf (".............\n");
    objcoe [i] = 0;
  }
}

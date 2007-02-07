//
// File: main.cpp
// Author: Pietro Belotti
// Purpose: test convexifier and function library
//
// (C) Pietro Belotti, 2006, Carnegie Mellon University
// This file is distributed under the Common Public License (CPL)
//

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#include <stdio.h>
#include <iostream>

#include <math.h>

//#include <expressions.h>
#include <CouenneTypes.h>
#include <CouenneProblem.h>
#include <CouenneCutGenerator.h>

#include <OsiSolverInterface.hpp>
#include <OsiSolverParameters.hpp>
#include <OsiClpSolverInterface.hpp>

// includes to read from a .nl file
#include <asl_pfgh.h>
#include "nlp2.h"
#include "getstub.h"
#include "r_opn.hd" /* for N_OPS */

// used within AMPL

#define CHR (char*)  

fint timing = 0;

static keyword keywds[] = { // must be alphabetical
   KW(CHR"timing", L_val, &timing, CHR"display timings for the run"),
};

static Option_Info Oinfo = { CHR"testampl", CHR"ANALYSIS TEST",
   CHR"concert_options", keywds, nkeywds, 0, CHR"ANALYSIS TEST" };

// create a problem from a file and add random convexification cuts

int main (int argc, char **argv) {

  // read .nl file from first argument //////////////////////////

  char *stub;

  // Create the ASL structure
  ASL_pfgh* asl = (ASL_pfgh*) ASL_alloc (ASL_read_pfgh);
  FILE *nl = NULL;
  stub = getstub   (&argv, &Oinfo);

  // Although very intuitive, we shall explain why the second argument
  // is passed with a minus sign: it is to tell the ASL to retrieve
  // the nonlinear information too.
  nl = jac0dim (stub, - (fint) strlen (stub));

  // Set options in the asl structure
  want_xpi0 = 1 | 2;  // allocate initial values for primal and dual if available
  obj_no = 0;         // always want to work with the first (and only?) objective

  // allocate space for initial values
  X0      = new real [n_var];
  havex0  = new char [n_var];
  pi0     = new real [n_con];
  havepi0 = new char [n_con];

  // read the rest of the nl file
  pfgh_read (nl, ASL_return_read_err | ASL_findgroups);


  // Example of use of a CouenneCutGenerator //////////////////////

  // create cut generator (the second parameter equal to true means
  // that only violated cuts should be generated)
  CouenneCutGenerator *cg = new CouenneCutGenerator 
    ((ASL_pfgh *) asl,  // ASL interface to the problem
     true,              // Do we want to get the violated cuts only?
     CURRENT_ONLY,      // how should convexification pick their sampling point 
     2);                // how many samples?

  // The constructor also sets all variable's bounds and initializes
  // them according to the information retrieved from the ASL_pfgh
  // structure. We can get it through the methods X(), Lb(), and
  // Ub(). Hence:

  // Allocate space for variables and bounds of convexification
  CouNumber *x = new CouNumber [cg -> getnvars ()];
  CouNumber *l = new CouNumber [cg -> getnvars ()];
  CouNumber *u = new CouNumber [cg -> getnvars ()];

  // get initial value (although in the beginning it is somewhat
  // useless)
  for (int i=0; i < cg -> getnvars (); i++) {
    x [i] = cg -> X  (i);
    l [i] = cg -> Lb (i);
    u [i] = cg -> Ub (i);
  }

  // print current state of the problem
  cg -> Problem () -> print (std::cout);

  // output initial values and bounds (bounds for original variables
  // are contained in the ASL structure. Initial value of original
  // variables can be provided by the ASL structure or is initialized
  // to 0. 
  //
  // The bounds (value) of each auxiliary variable is computed
  // as a function of the bounds (values) of the original variables 


  printf ("           x           lb             ub\n");
  for (int i=0; i < cg -> getnvars (); i++)
    printf ("%3d %12.3f %12.3f %12.3f\n", i, x [i], l [i], u [i]);


  // do whatever preprocessing on x, l, u within Bonmin 

  /* your code here */  

  // x,l,u should be relative to the original variables for the
  // first iteration only, and should be relative to all variables
  // starting from the second iteration on (to check violation).
  //  cg -> updateConv (x, l, u);

  // here we have obtained the first approximation. Solve the MILP
  // associated with this linearization. All cuts are available
  // through a call to "cg -> getCut (i)", which returns a pointer to
  // CouenneCut. Each of these objects are OsiRowCuts with an
  // additional CouNumber violation member, acc

  /*
  for (int i = 0; i < cg -> getncuts (); i++) {
    CouenneCut *cut = cg -> getCut (i);

    // ********************** your code here
  }  
  */

  // solving a convex problem with linearizations //////////////////////////////

  OsiClpSolverInterface *clp = new OsiClpSolverInterface;
  clp -> messageHandler () -> setLogLevel (1);

  OsiCuts cs;

  CglTreeInfo tree;

  cg -> generateCuts (*clp, cs, tree);

  int nvar0 = cg -> getnvars ();
  int ncon0 = cs . sizeRowCuts ();

  CoinBigIndex *start  = new CoinBigIndex [nvar0+1];
  int    *index  = (int    *) malloc (nvar0 * ncon0 * sizeof (int));
  double *coeff  = (double *) malloc (nvar0 * ncon0 * sizeof (double));
  double *objcoe = (double *) malloc (nvar0 * sizeof (double));
  double *rlb    = (double *) malloc (ncon0 * sizeof (double));
  double *rub    = (double *) malloc (ncon0 * sizeof (double));

  // fill start. First element is zero, for now just fill it with
  // number of nonzero elements for each column
  for (int i=nvar0; i--;)
    *start++ = 0;

  // empty the (nvar0+1)-th
  *start = 0;

  start -= nvar0;

  start++;

  // fill start with number of nonzero entries of each column
  for (int i = 0; i < ncon0; i++) {

    OsiRowCut *cut = cs. rowCutPtr (i);

    const int *ind   = cut -> row (). getIndices     ();
    int        nelem = cut -> row (). getNumElements (); 

    for (int j=0; j<nelem; j++)
      ++ (start [ind [j]]);
  }

  // cumulate elements of start []
  for (int i=1; i<nvar0; i++)
    start [i] += start [i-1];

  --start; // now elements in start are shifted right one position,
	   // and the contents are correct 

  for (int i=0; i < ncon0; i++) {

    OsiRowCut *cut = cs. rowCutPtr (i);

    const int    *ind   = cut -> row (). getIndices     ();
    const double *coe   = cut -> row (). getElements    ();
    int           nelem = cut -> row (). getNumElements (); 

    for (int j=0; j<nelem; j++) {

      int pos = start [ind [j]] ++;
      coeff [pos] = coe [j];
      index [pos] = i;
    }

    rlb [i] = cut -> lb ();
    rub [i] = cut -> ub ();
  }

  for (int i=nvar0; i; i--)
    start [i] = start [i-1];
  *start = 0;

  // set objective coefficients
  for (int i=0; i<nvar0; i++) 
    objcoe [i] = 0;

  objcoe [cg -> Problem () -> Obj (0) -> Body () -> Index ()] = 1;

  /*
  printf ("START: "); for (int i=0; i<=nvar0; i++) printf ("%2d ", start [i]); printf ("\n");
  printf ("Index: "); for (int i=0; i<start [nvar0]; i++) printf ("%2d ",index [i]); printf ("\n");
  printf ("coeff: "); for (int i=0; i<start [nvar0]; i++) printf ("%.3g ",coeff[i]); printf ("\n");
  printf ("objco: "); for (int i=0; i<nvar0; i++) printf ("%.3g ",  objcoe [i]); printf ("\n");
  printf ("lb:    "); for (int i=0; i<nvar0; i++)if(l[i]>-1e15) printf("%.3g ",l[i]);printf("\n");
  printf ("ub:    "); for (int i=0; i<nvar0; i++)if(u[i]< 1e15) printf("%.3g ",u [i]);printf ("\n");
  printf ("rlb:   "); for (int i=0; i<ncon0; i++)if(rlb[i]>-1e15)printf("%.3g ",rlb[i]);printf("\n");
  printf ("rub:   "); for (int i=0; i<ncon0; i++)if(rub[i]<1e15)printf("%.3g ",rub[i]);printf ("\n");
  */

  clp -> loadProblem (nvar0, ncon0, start, index, coeff, 
		      l, u, objcoe, rlb, rub);  
  delete [] start;

  free (index); free (coeff); free (objcoe);
  free (rlb);   free (rub);

  // ...and then call again the cut generation for a number of times:

  clp -> initialSolve ();
  printf ("-----------------------\n");

  // not more than 10 iterations of cut generation
  for (int round = 0; round < 5; round++) {

    for (int i = cs . sizeRowCuts (); i--;)
      cs. eraseRowCut (i);

    // x,l,u should be relative to the original variables for the
    // first iteration only, and should be relative to all variables
    // starting from the second iteration on (to check violation).
    //    cg -> updateConv (x, l, u);

    char fname [20];
    
    sprintf (fname, "round%d", round);

    //    clp -> writeMps (fname);
    clp -> writeLp  (fname);

    cg -> generateCuts (*clp, cs, tree);

    int ncuts = cs . sizeRowCuts ();

    /********************** some other code of yours here ***************************/
 
    // print generated cuts
    /*
    for (int i = 0; i < ncuts; i++)
      cs.rowCutPtr (i) -> print ();
    */

    // pseudo random sequence of bounds
    /*    
    for (int i=0; i < cg -> getnvars (); i++) {
      l [i] += 0.02;
      u [i] -= 0.4;
      x [i]  = 0.5 * (l [i] + u [i]);
    }
    */

    const OsiRowCut **newRowCuts = new const OsiRowCut * [ncuts];

    for (int i = 0; i < ncuts; i++)
      newRowCuts [i] = cs.rowCutPtr (i);
    /*
    for (int i = 0; i < ncuts; i++) {
      printf (" >> Adding cut ");
      newRowCuts [i] -> print ();
    }
    */
    clp -> applyRowCuts (ncuts, newRowCuts);

    delete [] newRowCuts;

    printf (">>> it. %4d: %4d new cuts. New obj.: %12.2f\n", 
	    round, ncuts, clp->getObjValue());

    clp -> resolve();  

    if (clp -> isAbandoned ())              
      {printf ("### ERROR: Numerical difficulties\n");break;}
    if (clp -> isProvenPrimalInfeasible ()) 
      {printf ("### WARNING: Problem infeasible\n");  break;}

    if (!ncuts) break;
  }

  delete clp;
  delete cg;

  delete [] x;  
  delete [] l;  
  delete [] u;

  delete [] X0;
  delete [] havex0;
  delete [] pi0;
  delete [] havepi0;

  ASL_free ((ASL **) &asl);

  return 0;
}

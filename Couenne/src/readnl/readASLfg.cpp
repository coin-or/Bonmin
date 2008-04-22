/*
 * Name:    readASLfg.cpp
 * Author:  Pietro Belotti
 * Purpose: re-read .nl file into a simpler ASL structure
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "asl.h"
#include "nlp.h"
#include "getstub.h"
#include "r_opn.hd" /* for N_OPS */

#define CHR (const char*)

fint timing = 0;

static keyword keywds[] = { // must be alphabetical
   KW(CHR"timing", L_val, &timing, CHR"display timings for the run"),
};

static Option_Info Oinfo = { CHR"testampl", CHR"ANALYSIS TEST",
   CHR"concert_options", keywds, nkeywds, 0, CHR"ANALYSIS TEST" };


ASL *readASLfg (char **argv) {

  // read .nl file from first argument //////////////////////////

  char *stub;

  // Create the ASL structure
  ASL* asl = (ASL*) ASL_alloc (ASL_read_fg);
  FILE *nl = NULL;
  stub = getstub (&argv, &Oinfo);

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
  fg_read (nl, ASL_return_read_err | ASL_findgroups);

  return asl;
}

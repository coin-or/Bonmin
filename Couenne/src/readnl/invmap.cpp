/*
 * Name: invmap.cpp
 * Author: Pietro Belotti
 * Purpose: create a bijection between ASL's efunc and integer to
 *          inversely map e->op fields into constant operators
 *
 * (C) 2006 Pietro Belotti, all rights reserved.
 * This code is distributed under the Common Public License.
 */

#include <stdlib.h>

#include <asl.h>
#include <opcode.hd>
#include <nlp.h>
#include <r_opn.hd>


/* couples an ASL function pointer with the relative operator constant */

typedef struct {
  efunc *fp;
  int    op;
} AslCouPair;


/* compare two AslCoupair's, used in qsort and bsearch below */

/* AW: 2007-06-11: changed b/c of problems with MSVC++ */
/* inline int pair_compare (const void *p1, const void *p2) { */
static int pair_compare (const void *p1, const void *p2) {

  /* FIX! weak cast for 64 bit machines */

  register long int f1 = (long int) (((AslCouPair *) p1) -> fp); 
  register long int f2 = (long int) (((AslCouPair *) p2) -> fp); 

  if      (f1 < f2) return -1;
  else if (f1 > f2) return  1;
  else return 0;
}


/* array of pairs (efunc2*, int) that relates all operators */

AslCouPair opmap [N_OPS];


/* binary search to get operator number from its efunc2* (the type of e->op) */

int getOperator (efunc *f) {

  static char first_call = 1;
  AslCouPair key, *res;

  /* FIX cast fo 64 bit machines */

  if (((long int) f <  N_OPS) && 
      ((long int) f > -N_OPS))
    return (long int) f;

  key.fp = f;

  if (first_call) { /* opmap is still empty, fill it using values from r_ops [] */

    register int i=0;
    register AslCouPair *ops = opmap;

    /* fill opmap vector with inverse correspondence pairs efunc -> int */
    while (i<N_OPS) {
      ops -> fp = r_ops [ops -> op = i++];
      ops++;
    }

    /* sort opmap for later use with bsearch */
    qsort (opmap, N_OPS, sizeof (AslCouPair), pair_compare);
    first_call = 0;
  }

  /* find int operator through binary search */
  res = (AslCouPair *) bsearch (&key, opmap, N_OPS, sizeof (AslCouPair), pair_compare);

  if (!res) 
    return -1;

  return res -> op;
}

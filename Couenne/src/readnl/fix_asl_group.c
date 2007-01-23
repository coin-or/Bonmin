/*
 * Name:    fix_asl_group.c
 * Author:  Pietro Belotti
 * Purpose: restore psg_groups to their original positions (i.e. inside expression tree)
 *
 * This file is licensed under the Common Public License (CPL)
 */


#include <asl_pfgh.h>
#include <nlp2.h>
#include <psinfo.h>
#include <opcode.hd>

#include <CouennePrecisions.h>


/* re-assign real expression to expr_n replacements of groups */

void fix_asl_group (psg_elem *g) {

  /*
   *  All psg_elem represent a sum of elements, and pfg*_read() reads
   *  the arguments of each unary operation as a group. That is,
   *  suppose we have the expression sin (exp (x0 + 3x1 + 3 +
   *  log(x1))) in a constraint or in an objective. The innermost
   *  expression, "x0 + 3x1 + 3 + log(x1)", is isolated from the tree
   *  and replaced with an expr_n, a constant expression whose value
   *  is updated at each call of AMPL's evaluation process.
   *
   *  The problem is that we need an expression tree and not a dead
   *  constant leaf, hence we have to restore the "x0 + 3x1 + 3 +
   *  log(x1)" where it was before.
   *
   *  Well, now that it is isolated, "x0 + 3x1 + 3 + log(x1)" is
   *  stored as a constant (g->g0 = 3), a set of g->nlin linear terms
   *  (here we have two, g->nlin=2) and a set of g->ns nonlinear terms
   *  (here g->ns=1). g->L is an array of g->nlin linpart's,
   *  containing an information on the coefficient (g->L[i].fac) and
   *  index (g -> L[j]. v.i) of the term. Finally, g->E is an array of
   *  expression, here containing only one cell with the expression
   *  tree for log(x1). 
   *
   *  We need to create an OPSUMLIST expression whose terms are all
   *  the isolated terms in "x0 + 3x1 + 3 + log(x1)".
   *
   */

  int j, nterms = g->ns + g->nlin;

  expr **arglist;
  expr  *sum;

  if (fabs (g -> g0) > COUENNE_EPS) 
    nterms++;

  /* argument list and expression of OPSUMLIST */
  arglist = (expr **) malloc (nterms * sizeof (expr *));
  sum     = (expr  *) malloc (sizeof (expr));

  sum -> L.ep = arglist;
  sum -> R.ep = arglist + nterms;

  /* add constant term only if nonzero */
  if (fabs (g -> g0) > COUENNE_EPS) {
    expr_n *const_term = (expr_n *)
      (arglist [--nterms] = (expr *) malloc (sizeof (expr_n)));

    const_term -> op = (efunc_n *)(long) OPNUM;
    const_term -> v = g -> g0;
  }

  /* add nonlinear terms */
  for (j = 0; j < g->ns; j++)
    arglist [--nterms] = g -> E[j].D.e;

  /* now the linear part */
  for (j = 0; j < g->nlin; j++) {

    /* coefficient */
    real coef = g -> L[j]. fac;
    expr *lin = (expr *) malloc (sizeof (expr));

    if (fabs (coef - 1.0) > COUENNE_EPS) { 

      /* expression is ax, a!=1 */

      expr_n *num = (expr_n *) malloc (sizeof (expr_n));
      num -> op = (efunc_n *)(long) OPNUM;
      num -> v  = coef;

      lin -> L.e = (expr *) num;

      lin -> R.e = (expr *) malloc (sizeof (expr));
      lin -> R.e -> op = (efunc2 *)(long) OPVARVAL;
      lin -> R.e -> a = g -> L[j]. v.i;

      lin -> op = (efunc2 *)(long) OPMULT;

    } else {

      /* expression is x */

      /*      lin = (expr *) malloc (sizeof (expr));*/
      lin -> op = (efunc2 *)(long) OPVARVAL;
      lin -> a = g -> L[j]. v.i;
    }

    arglist [--nterms] = lin;
  }

  sum -> op = (efunc2 *)(long) OPSUMLIST;
  g -> ge -> L.e = sum;
}


/* free restored group (since not taken care of within ASL_free) */

void free_asl_group (psg_elem *g) {

  int j;

  /* free constant expr_n part */
  if (fabs (g -> g0) > COUENNE_EPS)
    free (g -> ge -> L.e -> L.ep [g->ns + g->nlin]);

  /* free linear parts */
  for (j = 0; j < g->nlin; j++) {
    if (fabs (g -> L[j]. fac - 1.0) > COUENNE_EPS) {
      free (g -> ge -> L.e -> L.ep [j] -> L.e);
      free (g -> ge -> L.e -> L.ep [j] -> R.e);
    }
    free (g -> ge -> L.e -> L.ep [j]);
  }

  /* free the whole OPSUMLIST */

  free (g -> ge -> L.e -> L.ep);
  /*  free (g -> ge -> L.e); */
  g -> ge -> L.e = (expr*)&g -> esum;
  /*  g -> ge -> L.e = g -> ge.esum; */
}

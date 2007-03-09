/*
 * Name:    readnl.cpp
 * Author:  Pietro Belotti
 * Purpose: define a reader for .nl files. Adapted from ampl2ev3 by L. Liberti and S. Galli 
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneProblem.h>

#include <CouenneTypes.h>

#include <exprSum.h>
#include <exprMul.h>
#include <exprClone.h>
#include <exprGroup.h>

#include "asl.h"
#include "nlp.h"
#include "getstub.h"
#include "opcode.hd"

#define OBJ_DE    ((const ASL_fg *) asl) -> I.obj_de_
#define CON_DE    ((const ASL_fg *) asl) -> I.con_de_
#define OBJ_sense ((const ASL_fg *) asl) -> i.objtype_


// check if an expression is a null pointer or equals zero

inline bool is_expr_zero (expr* e)
  {return ((e==NULL) || (((((long int) e->op) == OPNUM) && 
			  (fabs (((expr_n *)e) -> v)  < COUENNE_EPS) 
			  //  && (fabs (e -> dL) < COUENNE_EPS)
			  // *** CHECK THIS! dL is the derivative
			  )));} 


// Reads a MINLP from an AMPL .nl file through the ASL methods

int CouenneProblem::readnl (const ASL *asl) {

  int n_intvar = niv + nbv + nlvbi + nlvci + nlvoi;

  // create continuous variables
  for (int i = n_var-n_intvar; i--;)
    addVariable (false);

  // create integer variables
  for (int i = n_intvar; i--;)
    addVariable (true);

  //  printf ("readln: added %d (%d + %d) variables\n", n_var, n_var - n_intvar, n_intvar);

  // create room for problem's variables and bounds
  x_  = (CouNumber *) malloc (n_var * sizeof (CouNumber));
  lb_ = (CouNumber *) malloc (n_var * sizeof (CouNumber));
  ub_ = (CouNumber *) malloc (n_var * sizeof (CouNumber));


  // objective functions /////////////////////////////////////////////////////////////

  for (int i = 0; i < n_obj; i++) {

    ////////////////////////////////////////////////
    int nterms = 0;

    // count nonzero terms in linear part
 
    for (ograd *objgrad = Ograd [i];
	 objgrad;
	 objgrad = objgrad -> next)
      if (fabs (objgrad -> coef) > COUENNE_EPS)
	nterms++;

    expression *body;

    expression *nl = nl2e (OBJ_DE [i] . e);

    if (nterms) {

      int       *index = new int       [nterms+1];
      CouNumber *coeff = new CouNumber [nterms];

      for (ograd *objgrad = Ograd [i]; objgrad; objgrad = objgrad -> next)
	if (fabs (objgrad -> coef) > COUENNE_EPS) {

	  *index++ = objgrad -> varno;
	  *coeff++ = objgrad -> coef;
	}

      *index = -1;

      index -= nterms;
      coeff -= nterms;

      expression **nll = new expression * [1];

      *nll = nl;

      body = new exprGroup (objconst (i), index, coeff, nll, 1);

      delete [] index;
      delete [] coeff;
    } else 
      if (fabs (objconst (i) > COUENNE_EPS))
	body = new exprSum (nl, new exprConst (objconst (i)));
      else body = nl;

    ///////////////////////////////////////////////////

    expression *subst = body -> simplify ();
    if (subst) body = subst;

    // ThirdParty/ASL/solvers/asl.h, line 336: 0 is minimization, 1 is maximization
    addObjective (body, (OBJ_sense [i] == 0) ? "min" : "max");
  }


  // constraints ///////////////////////////////////////////////////////////////////

  int          *nterms = new int           [n_con];

  // allocate space for argument list of all constraints' summations
  // of linear and nonlinear terms

  // init array with # terms of each constraint
  for (int i = n_con; i--;) 
    *nterms++ = 0;
  nterms -= n_con;

  cgrad *congrad;

  // count all linear terms
  if (A_colstarts && A_vals)
    for (register int j = A_colstarts [n_var]; j--;) {

      real coeff = A_vals [j];

      if (fabs (coeff) > COUENNE_EPS)
	nterms [A_rownos [j]] ++;
    }
  else {                             // Constraints' linear info is stored in Cgrad
    for (register int i = 0; i < n_con; i++)
      for (congrad = Cgrad [i]; 
	   congrad; 
	   congrad = congrad -> next) 
	if (fabs (congrad -> coef) > COUENNE_EPS) 
	  nterms [i] ++;
  }


  // vectors of the linear part
  CouNumber **coeff = new CouNumber * [n_con];
  int       **index = new int       * [n_con];

  for (register int i = n_con; i--;) 
    *index++ = NULL;

  index -= n_con;


  // set linear terms

  if (A_colstarts && A_vals)         // Constraints' linear info is stored in A_vals
    for (int j = 0; j < n_var; j++)
      for (register int i = A_colstarts [j], k = A_colstarts [j+1] - i; k--; i++) {

	int rowno = A_rownos [i],
	    nt    = nterms [rowno] --;

	CouNumber **cline = coeff + rowno;
	int       **iline = index + rowno;

	if (*iline==NULL) {
	  *cline = new CouNumber [nt];
	  *iline = new int       [nt+1];
	  (*iline) [nt] = -1;
	}

	(*cline) [--nt] = A_vals [i];
	(*iline)   [nt] = j;

      }
  else {                             // Constraints' linear info is stored in Cgrad
    for (int i=0; i < n_con; i++) {

      int nt = nterms [i];

      CouNumber **cline = coeff + i;
      int       **iline = index + i;

      *cline = new CouNumber [nt];
      *iline = new int       [nt+1];
      (*iline) [nt] = -1;

      for (congrad = Cgrad [i]; congrad; congrad = congrad -> next) 
	if (fabs (congrad -> coef) > COUENNE_EPS) {
	  (*cline) [--nt] = congrad -> coef;
	  (*iline)   [nt] = congrad -> varno;
	}
    }
  }

  // set constraints' bound and sign and store nonlinear part ///////////////////////////////

  for (int i = 0; i < n_con; i++) {

    enum con_sign sign;
    double lb, ub;

    if (Urhsx) {
      lb = LUrhs [i];
      ub = Urhsx [i];
    } else {
      int j = 2*i;
      lb = LUrhs [j];
      ub = LUrhs [j+1];
    }

    // set constraint sign
    if (lb > negInfinity + 1) 
      if (ub < Infinity - 1) 
	sign = COUENNE_RNG;
      else sign = COUENNE_GE;
    else sign = COUENNE_LE;

    // this is an equality constraint  
    if (fabs (lb - ub) < COUENNE_EPS)
      sign = COUENNE_EQ;

    expression *body;

    expression **nll = new expression * [1];
    *nll = nl2e (CON_DE [i] . e);

    if (index [i]) 
      body = new exprGroup (0, index [i], coeff [i], nll, 1);
    else {
      body = *nll;
      delete [] nll;
    }

    expression *subst = body -> simplify ();
    if (subst) body = subst;

    // add them (and set lower-upper bound)
    switch (sign) {
	
    case COUENNE_EQ:  addEQConstraint  (body, new exprConst (ub)); break;
    case COUENNE_LE:  addLEConstraint  (body, new exprConst (ub)); break;
    case COUENNE_GE:  addGEConstraint  (body, new exprConst (lb)); break;
    case COUENNE_RNG: addRNGConstraint (body, new exprConst (lb), 
					      new exprConst (ub)); break;
    default: printf ("could not recognize constraint\n"); return -1;
    }

    delete [] index [i];
    delete [] coeff [i];
  }

  delete [] index;
  delete [] coeff;

  // lower and upper bounds ///////////////////////////////////////////////////////////////

  if (LUv) {

    real *Uvx_copy = Uvx;

    if (!Uvx_copy)
      for (register int i=0; i<n_var; i++) {

	register int j = 2*i;

	lb_ [i] = LUv [j];
	ub_ [i] = LUv [j+1];
      }
    else
      for (register int i=n_var; i--;) {
	lb_ [i] = LUv [i];
	ub_ [i] = Uvx_copy [i];
      }

  } else
    for (register int i=n_var; i--;) {
      lb_ [i] = - COUENNE_INFINITY;
      ub_ [i] =   COUENNE_INFINITY;
    }

  // initial x ////////////////////////////////////////////////////////////////////

  for (register int i=n_var; i--;) 
    if (X0 && havex0 [i])
      x_ [i] = X0 [i]; 
    else {

      CouNumber x, l = lb_ [i], u = ub_ [i];

      if      (l < - COUENNE_INFINITY + 1)
	if    (u >   COUENNE_INFINITY - 1)  x = 0;
	else                                x = u;
      else if (u >   COUENNE_INFINITY - 1)  x = l;
      else                                  x = 0.5 * (l+u);

      x_ [i] = x;
    }

  // update expression internal data

  expression::update (x_, lb_, ub_);

  delete [] nterms;

  return 0;
}

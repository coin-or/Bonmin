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

#include "asl_pfgh.h"
#include "nlp2.h"
#include "getstub.h"
#include "opcode.hd"

#define OBJ_DE ((const ASL_pfgh *) asl) -> I.obj2_de_
#define CON_DE ((const ASL_pfgh *) asl) -> I.con2_de_
#define OBJ_sense ((const ASL_pfgh *) asl) -> i.objtype_


// check if an expression is a null pointer or equals zero

inline bool is_expr_zero (expr* e)
  {return ((e==NULL) || (((((long int) e->op) == OPNUM) && 
			  (fabs (((expr_n *)e) -> v)  < COUENNE_EPS) 
			  //  && (fabs (e -> dL) < COUENNE_EPS)
			  // *** CHECK THIS! dL is the derivative
			  )));} 


// replaces group in original position within expression tree

extern "C" {
  void fix_asl_group  (psg_elem *);
  void free_asl_group (psg_elem *);
}


// Reads a MINLP from an AMPL .nl file through the ASL methods

int CouenneProblem::readnl (const ASL_pfgh *asl) {

  int n_intvar = niv + nbv + nlvbi + nlvci + nlvoi;

  // create discrete variables
  for (int i = n_var-n_intvar; i--;)
    addVariable (false);

  // create continuous variables
  for (int i = n_intvar; i--;)
    addVariable (true);

  // create room for problem's variables and bounds
  x_  = (CouNumber *) malloc (n_var * sizeof (CouNumber));
  lb_ = (CouNumber *) malloc (n_var * sizeof (CouNumber));
  ub_ = (CouNumber *) malloc (n_var * sizeof (CouNumber));


  // fix groups //////////////////////////////////////////////////////////////////////

  // objective functions' groups

  for (int k=0; k<n_obj; k++)
    if ((asl -> P. ops [k]. nb) || (asl -> P. ops [k]. ng)) {

      int ngroups = asl -> P. ops [k]. ng;
      for (int i=0; i<ngroups; i++)
	fix_asl_group (&(asl -> P. ops [k]. g [i]));
    }

  // constraints' groups

  for (int k=0; k<n_con; k++)
    if ((asl -> P. cps [k]. nb) || (asl -> P. cps [k]. ng)) {

      int ngroups = asl -> P. cps [k]. ng;
      for (int i=0; i<ngroups; i++)
	fix_asl_group (&(asl -> P. cps [k]. g [i]));
    }

  // objective functions /////////////////////////////////////////////////////////////

  for (int i = 0; i < n_obj; i++) {

    int nterms = 0;

    real z0 = objconst (i);

    // count nonzero terms in objective function (constant, linear
    // and nonlinear part)
 
    if (fabs (z0) > COUENNE_EPS) nterms++;

    for (ograd *objgrad = Ograd [i]; 
	 objgrad; 
	 objgrad = objgrad -> next) 
      if (fabs (objgrad -> coef) > COUENNE_EPS) 
	nterms++;

    if ((OBJ_DE [i] . e) && (!is_expr_zero (OBJ_DE [i] . e)))
      nterms++;

    if (!nterms) // strange, but no (linear or nonlinear) terms in here 
      continue;

    // now we can fill in the objective function terms

    expression **al = new expression * [nterms];

    int term = 0;

    if (fabs (z0) > COUENNE_EPS) 
      al [term++] = new exprConst (z0);

    ograd *objgrad;

    for (objgrad = Ograd [i]; 
	 objgrad; 
	 objgrad = objgrad -> next) 
      if (fabs (objgrad -> coef) > COUENNE_EPS) 
	al [term++] = new exprMul (new exprConst (objgrad -> coef),
				   new exprClone (variables_ [objgrad -> varno]));

    if (   (OBJ_DE [i] . e)
	   && (!is_expr_zero (OBJ_DE [i] . e)))
      al [term++] = nl2e (OBJ_DE [i] . e);

    expression *body;

    // if sum only has one term, this objective is rather that term than a sum...
    if (term==1) {
      body = al [0]; 
      delete [] al;
    }
    else 
      body = new exprSum (al, term); 

    expression *subst = body -> simplify ();
    if (subst) body = subst;

    // ThirdParty/ASL/solvers/asl.h, line 336: 0 is minimization, 1 is maximization
    addObjective (body, (OBJ_sense [i] == 0) ? "min" : "max");
  }


  // constraints ///////////////////////////////////////////////////////////////////

  expression ***alists = new expression ** [n_con];
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
    for (int j = A_colstarts [n_var]; j--;) {

      real coeff = A_vals [j];

      if (fabs (coeff) > COUENNE_EPS)
	nterms [A_rownos [j]] ++;
    }
  else {                             // Constraints' linear info is stored in Cgrad
    for (int i = 0; i < n_con; i++)
      for (congrad = Cgrad [i]; 
	   congrad; 
	   congrad = congrad -> next) 
	if (fabs (congrad -> coef) > COUENNE_EPS) 
	  nterms [i] ++;
  }


  // add possible nonlinear term (at most one)
  for (int i = 0; i < n_con; i++) 
    if (!is_expr_zero (CON_DE [i] . e))
      nterms [i] ++;

  // reserve space for argument lists of sums
  for (int i = 0; i < n_con; i++)
    alists [i] = new expression * [nterms [i]];

  // set linear terms

  if (A_colstarts && A_vals)         // Constraints' linear info is stored in A_vals
    for (int j = 0; j < n_var; j++)
      for (int i = A_colstarts [j]; i < A_colstarts [j + 1]; i++) {

	real coeff = A_vals[i];

	if (fabs (coeff) > COUENNE_EPS)
	  if (fabs (coeff - 1.0) > COUENNE_EPS)
	    *(alists [A_rownos [i]])++ = 
	      new exprMul (new exprConst (coeff), new exprClone (variables_ [j]));
	  else
	    *(alists [A_rownos [i]])++ = 
	      new exprClone (variables_ [j]);
      }
  else {                             // Constraints' linear info is stored in Cgrad
    for (int i = 0; i < n_con; i++)
      for (congrad = Cgrad [i]; 
	   congrad; 
	   congrad = congrad -> next) 
	if (fabs (congrad -> coef) > COUENNE_EPS) 
	  *(alists [i])++ = new exprMul (new exprConst (congrad -> coef),
					 new exprClone (variables_ [congrad -> varno]));
  }


  // set constraints' bound, sign and store nonlinear part ///////////////////////////////

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
  
    if (fabs (lb - ub) < COUENNE_EPS)
      sign = COUENNE_EQ;

    // set nonlinear part
    if (!is_expr_zero (CON_DE [i] . e))
      *(alists [i])++ = nl2e (CON_DE [i] . e);

    // re-align alists [i] to original value
    alists [i] -= nterms [i];

    if (!(nterms [i])) 
      delete [] alists [i];
    else {

      expression *body;

      // create sum with alists [i] as argument list
      if (nterms [i] == 1) body = alists [i] [0];
      else                 body = new exprSum (alists [i], nterms [i]);	

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
    }
  }

  // lower and upper bounds ///////////////////////////////////////////////////////////////

  if (LUv) {

    real *Uvx_copy = Uvx;

    if (!Uvx_copy)
      for (int i=0; i<n_var; i++) {

	int j = 2*i;

	lb_ [i] = LUv [j];
	ub_ [i] = LUv [j+1];
      }
    else
      for (int i=n_var; i--;) {
	lb_ [i] = LUv [i];
	ub_ [i] = Uvx_copy [i];
      }

  } else
    for (int i=n_var; i--;) {
      lb_ [i] = - COUENNE_INFINITY;
      ub_ [i] =   COUENNE_INFINITY;
    }

  // initial x ////////////////////////////////////////////////////////////////////

  for (int i=n_var; i--;) 
    if (X0 && havex0 [i])
      x_ [i] = X0 [i]; 
    else {
      //    for (int i=n_var; i--;) {

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
  delete [] alists;

  // restore groups for use by Ipopt: constraints

  for (int k=0; k<n_con; k++)
    if ((asl -> P. cps [k]. nb) || (asl -> P. cps [k]. ng)) {

      int ngroups = asl -> P. cps [k]. ng;
      for (int i=0; i<ngroups; i++)
	free_asl_group (asl -> P. cps [k]. g + i);
    }

  // and objectives

  for (int k=0; k<n_obj; k++)
    if ((asl -> P. ops [k]. nb) || (asl -> P. ops [k]. ng)) {

      int ngroups = asl -> P. ops [k]. ng;
      for (int i=0; i<ngroups; i++)
	free_asl_group (asl -> P. ops [k]. g + i);
    }

  return 0;
}

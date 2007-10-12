/*
 * Name:    readnl.cpp
 * Author:  Pietro Belotti
 * Purpose: define a reader for .nl files. Adapted from ampl2ev3 by L. Liberti and S. Galli 
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneProblem.hpp>

#include <CouenneTypes.hpp>

#include <exprSum.hpp>
#include <exprMul.hpp>
#include <exprClone.hpp>
#include <exprGroup.hpp>

#include "asl.h"
#include "nlp.h"
#include "getstub.h"
#include "opcode.hd"

#define OBJ_DE    ((const ASL_fg *) asl) -> I.obj_de_
#define VAR_E     ((const ASL_fg *) asl) -> I.var_e_
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

  // # discrete variables
  int n_intvar = niv + nbv + nlvbi + nlvci + nlvoi;

  // number of defined variables (aka common expressions)
  ndefined_ = como + comc + comb + como1 + comc1; 

  // create continuous variables
  for (int i = n_var - n_intvar; i--;)
    addVariable (false);

  // create integer variables
  for (int i = n_intvar; i--;)
    addVariable (true);

  nOrig_ = n_var;

  // create expression set for binary search
  auxSet_ = new std::set <exprAux *, compExpr>;

  // create room for problem's variables and bounds
  x_  = (CouNumber *) malloc ((n_var + ndefined_) * sizeof (CouNumber));
  lb_ = (CouNumber *) malloc ((n_var + ndefined_) * sizeof (CouNumber));
  ub_ = (CouNumber *) malloc ((n_var + ndefined_) * sizeof (CouNumber));

  for (int i = n_var + ndefined_; i--;) {
    x_  [i] =  0;
    lb_ [i] = -COUENNE_INFINITY;
    ub_ [i] =  COUENNE_INFINITY;
  }

  // common expressions (or defined variables) ///////////////////////////////////////

  /*printf ("c_vars_ = %d\n", ((const ASL_fg *) asl) -> i.c_vars_ );
  printf ("comb_ = %d\n",   ((const ASL_fg *) asl) -> i.comb_  );
  printf ("combc_ = %d\n",  ((const ASL_fg *) asl) -> i.combc_ );
  printf ("comc1_ = %d\n",  ((const ASL_fg *) asl) -> i.comc1_ );
  printf ("comc_ = %d\n",   ((const ASL_fg *) asl) -> i.comc_  );
  printf ("como1_ = %d\n",  ((const ASL_fg *) asl) -> i.como1_ );
  printf ("como_ = %d\n",   ((const ASL_fg *) asl) -> i.como_  );*/

  // Each has a linear and a nonlinear part (credit to Dominique
  // Orban: http://www.gerad.ca/~orban/drampl/def-vars.html)

  for (int i = 0; i < como + comc + comb; i++) {

    struct cexp *common = ((const ASL_fg *) asl) -> I.cexps_ + i;
    expression *nle = nl2e (common -> e);

    //    printf ("cexp  %d: ", i); nle -> print ();  printf (" ||| ");

    int nlin = common -> nlin;  // Number of linear terms

    if (nlin > 0) {

      int       *index = new int       [nlin+1];
      CouNumber *coeff = new CouNumber [nlin];

      linpart *L = common -> L;

      for (int j = 0; j < nlin; j++ ) {
	//	vp = (expr_v *)((char *)L->v.rp - ((char *)&ev.v - (char *)&ev));
	//	Printf( " %-g x[%-d]", L->fac, (int)(vp - VAR_E) );	
	coeff [j] = L [j]. fac;
	index [j] = (expr_v *) (L [j].v.rp) - VAR_E;
	/*	Printf( " %-g x[%-d]", L [j]. fac, 
		(expr_v *) (L [j].v.rp) - VAR_E //((const ASL_fg *) asl) -> I.cexps_
		//L [j]. v.i
		);*/
      }

      index [nlin] = -1;

      expression **al = new expression * [1];
      *al = nle;

      exprGroup *eg = new exprGroup (0, index, coeff, al, 1);
      commonexprs_ . push_back (eg);
    } 
    else commonexprs_ . push_back (nle);
    //    printf ("\n");
  }

  for (int i = 0; i < como1 + comc1; i++) {

    struct cexp1 *common = ((const ASL_fg *) asl) -> I.cexps1_ + i;
    expression *nle = nl2e (common -> e);

    //    printf ("cexp1 %d: ", i); nle -> print ();  printf (" ||| ");

    int nlin = common -> nlin;  // Number of linear terms

    if (nlin > 0) {

      int       *index = new int       [nlin+1];
      CouNumber *coeff = new CouNumber [nlin];

      linpart *L = common -> L;

      for (int j = 0; j < nlin; j++ ) {
	//	vp = (expr_v *)((char *)L->v.rp - ((char *)&ev.v - (char *)&ev));
	coeff [j] = L [j]. fac;
	index [j] = (expr_v *) (L [j].v.rp) - VAR_E;
	/*	Printf( " %-g x[%-d]", L [j]. fac, 
		(expr_v *) (L [j].v.rp) - VAR_E //((const ASL_fg *) asl) -> I.cexps_
		//L [j]. v.i
		);*/
      }

      index [nlin] = -1;

      expression **al = new expression * [1];
      *al = nle;

      exprGroup *eg = new exprGroup (0, index, coeff, al, 1);
      commonexprs_ . push_back (eg);
    } 
    else commonexprs_ . push_back (nle);
    //    printf ("\n");
    //    addAuxiliary (nl2e (((const ASL_fg *) asl) -> I.cexps1_ [i] . e));
  }

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

      if (nl -> code () == COU_EXPRSUM)
	body = new exprGroup (0, index, coeff, nl -> ArgList (), nl -> nArgs ());
      else {

	expression **nll = new expression * [1];

	*nll = nl;

	//body = new exprGroup (objconst (i), index, coeff, nll, 1);

	// apparently, objconst (i) is included in the obj expression
	body = new exprGroup (0, index, coeff, nll, 1);
      }

      delete [] index;
      delete [] coeff;

    } else
      // apparently, objconst (i) is included in the obj expression
      body = nl;
      //if (fabs (objconst (i) > COUENNE_EPS))
      //body = new exprSum (nl, new exprConst (objconst (i)));
      //else 

    ///////////////////////////////////////////////////

    expression *subst = body -> simplify ();
    if (subst) {
      delete body; // VALGRIND
      body = subst;
    }

    // ThirdParty/ASL/solvers/asl.h, line 336: 0 is minimization, 1 is maximization
    addObjective (body, (OBJ_sense [i] == 0) ? "min" : "max");
  }


  // constraints ///////////////////////////////////////////////////////////////////

  int *nterms = new int [n_con];

  // allocate space for argument list of all constraints' summations
  // of linear and nonlinear terms

  // init array with # terms of each constraint
  for (int i = n_con; i--;) 
    *nterms++ = 0;
  nterms -= n_con;

  cgrad *congrad;

  // count all linear terms
  if (A_colstarts && A_vals)         // Constraints' linear info is stored in A_vals
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
    if (lb > negInfinity)
      if (ub < Infinity) sign = COUENNE_RNG;
      else               sign = COUENNE_GE;
    else                 sign = COUENNE_LE;

    // this is an equality constraint  
    if (fabs (lb - ub) < COUENNE_EPS)
      sign = COUENNE_EQ;

    expression *body;

    expression **nll = new expression * [1];
    *nll = nl2e (CON_DE [i] . e);

    if (index [i] && (*(index [i]) >= 0)) 
      if ((*nll) -> code () == COU_EXPRSUM)
	body    = new exprGroup (0., index [i], coeff [i], (*nll) -> ArgList (), (*nll) -> nArgs ());
      else body = new exprGroup (0., index [i], coeff [i], nll, 1);
    else {
      body = *nll;
      delete [] nll;
    }

    expression *subst = body -> simplify ();
    if (subst) {
      delete body; // VALGRIND
      body = subst;
    }

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

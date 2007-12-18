/*
 * Name:    nl2e.cpp
 * Author:  Pietro Belotti
 * Purpose: converts a nl expression into a Couenne expression
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouenneTypes.hpp"

#include "exprAbs.hpp"
#include "exprSum.hpp"
#include "exprSub.hpp"
#include "exprMul.hpp"
#include "exprDiv.hpp"
#include "exprInv.hpp"
#include "exprSin.hpp"
#include "exprPow.hpp"
#include "exprClone.hpp"
#include "exprLog.hpp"
#include "exprOpp.hpp"
#include "exprCos.hpp"
#include "exprExp.hpp"

#include "asl.h"
#include "nlp.h"
#include "opcode.hd"


// get ASL op. code relative to function pointer passed as parameter 
int getOperator (efunc *);


// warning for non-implemented functions -- return 0 constant expression
//expression *notimpl (const std::string &fname) {
void notimpl (const std::string &fname) {
  std::cerr << "*** Error: " << fname << " not implemented" << std::endl;
  exit (-1);
}


// converts an AMPL expression (sub)tree into an expression* (sub)tree

expression *CouenneProblem::nl2e (expr *e) {

  switch (getOperator (e -> op)) {

  case OPPLUS:  return new exprSum (nl2e (e -> L.e), nl2e (e -> R.e));
  case OPMINUS: return new exprSub (nl2e (e -> L.e), nl2e (e -> R.e));
  case OPMULT:  return new exprMul (nl2e (e -> L.e), nl2e (e -> R.e));
  case OPDIV:   return new exprDiv (nl2e (e -> L.e), nl2e (e -> R.e));
  case OPREM:   notimpl ("remainder");
  case OPPOW:   return new exprPow (nl2e (e -> L.e), nl2e (e -> R.e));
  case OPLESS:  notimpl ("less");
  case MINLIST: notimpl ("min");
  case MAXLIST: notimpl ("max");
  case FLOOR:   notimpl ("floor");
  case CEIL:    notimpl ("ceil");
  case ABS:     return new exprAbs (nl2e (e -> L.e));
  case OPUMINUS:return new exprOpp (nl2e (e -> L.e));
    //  return new exprOpp (nl2e (e -> L.e -> L.e));
  case OPIFnl:  { notimpl ("ifnl");

  // see ASL/solvers/rops.c, IfNL
  }

  case OP_tanh: notimpl ("tanh");
  case OP_tan: {
    expression *arg;
    arg = nl2e (e -> L.e);
    return new exprDiv (new exprSin (arg), new exprCos (new exprClone (arg)));
  }
  case OP_sqrt:    return new exprPow (nl2e (e -> L.e), new exprConst (0.5));
  case OP_sinh:    return new exprMul (new exprConst (0.5),
				       new exprSub (new exprExp (nl2e (e -> L.e)),
						    new exprExp (new exprOpp (nl2e (e -> L.e)))));
  case OP_sin:     return new exprSin (nl2e (e -> L.e));
  case OP_log10:   return new exprMul (new exprConst (1.0 / log (10.0)), 
				       new exprLog (nl2e (e -> L.e)));
  case OP_log:     return new exprLog (nl2e (e -> L.e));
  case OP_exp:     return new exprExp (nl2e (e -> L.e));
  case OP_cosh:    return new exprMul (new exprConst (0.5),
				       new exprSum (new exprExp (nl2e (e -> L.e)),
						    new exprExp (new exprOpp (nl2e (e -> L.e)))));

  case OP_cos:   return new exprCos (nl2e (e -> L.e));
  case OP_atanh: notimpl ("atanh");
  case OP_atan2: notimpl ("atan2");
  case OP_atan:  notimpl ("atan");
  case OP_asinh: notimpl ("asinh");
  case OP_asin:  notimpl ("asin");
  case OP_acosh: notimpl ("acosh");
  case OP_acos:  notimpl ("acos");

  case OPSUMLIST: {
    register int i=0;
    expression **al = new expression * [(e->R.ep - e->L.ep)];
    for (expr **ep = e->L.ep; ep < e->R.ep; ep++)
      al [i++] = nl2e (*ep);
    return new exprSum (al, i);
  }
  case OPintDIV: notimpl ("intdiv");
  case OPprecision: notimpl ("precision");
  case OPround:  notimpl ("round");
  case OPtrunc:  notimpl ("trunc");

  case OP1POW: return new exprPow (nl2e (e -> L.e), 
				   new exprConst (((expr_n *)e->R.e)->v));
  case OP2POW: return new exprPow (nl2e (e -> L.e), 
				   new exprConst (2.));
  case OPCPOW: return new exprPow (new exprConst (((expr_n *)e->L.e)->v),
				   nl2e (e -> R.e));
  case OPFUNCALL: notimpl ("function call");
  case OPNUM:     return new exprConst (((expr_n *)e)->v);
  case OPPLTERM:  notimpl ("plterm");
  case OPIFSYM:   notimpl ("ifsym");
  case OPHOL:     notimpl ("hol");
  case OPVARVAL:  {

    int j = ((expr_v *) e) -> a, 
        d = nVars () - j;

    if (j >= nOrig_) j--; // CHECK! 

    d = nVars () - j;

    // is index above number of variables?
    //if (j >= nvars)

    while (d++ <= 0)
      addVariable (false);

    //printf ("indexD = %d\n", j);
    //printf ("Couenne, warning: unknown variable x%d (>%d+%d=%d), returning new variable.\n",
    //	j, nvars, nAuxs (), nvars + nAuxs ());
    //	exit (-1);
    // TODO: aux_ [...] may not be filled yet

    return new exprClone (variables_ [j]);
  }

  default: 
    printf ("ERROR: unknown operator (address %x), aborting.\n", (long int) e -> op); 
    exit (-1);
    //return new exprConst (0);
  }

  return new exprConst (0.);
}

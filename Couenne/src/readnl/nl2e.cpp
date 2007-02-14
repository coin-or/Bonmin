/*
 * Name:    nl2e.cpp
 * Author:  Pietro Belotti
 * Purpose: converts a nl expression into a Feline expression
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>

#include <exprAbs.h>
#include <exprSum.h>
#include <exprSub.h>
#include <exprMul.h>
#include <exprDiv.h>
#include <exprInv.h>
#include <exprSin.h>
#include <exprPow.h>
#include <exprClone.h>
#include <exprLog.h>
#include <exprOpp.h>
#include <exprCos.h>
#include <exprExp.h>

#include <asl_pfgh.h>
#include <nlp2.h>
#include <opcode.hd>


// get ASL op. code relative to function pointer passed as parameter 

extern "C" {
  int getOperator (efunc2 *);
}


// converts an AMPL expression (sub)tree into an expression* (sub)tree

expression *CouenneProblem::nl2e (expr2 *e) {

  expression **al;
  expr **ep;
  expression *arg;

  switch (getOperator (e -> op)) {

  case OPPLUS:   return new exprSum (nl2e (e -> L.e), nl2e (e -> R.e));
  case OPMINUS:  return new exprSub (nl2e (e -> L.e), nl2e (e -> R.e -> L.e));
  case OPMULT:   return new exprMul (nl2e (e -> L.e), nl2e (e -> R.e));
  case OPDIV:    return new exprDiv (nl2e (e -> L.e), nl2e (e -> R.e));
  case OPREM:   printf ("remainder not implemented\n"); return new exprConst (0);
  case OPPOW:    return new exprPow (nl2e (e -> L.e), nl2e (e -> R.e));
  case OPLESS:  printf ("less not implemented\n"); return new exprConst (0);
  case MINLIST: printf ("min not implemented\n"); return new exprConst (0);
  case MAXLIST: printf ("max not implemented\n"); return new exprConst (0);
  case FLOOR:   printf ("floor not implemented\n"); return new exprConst (0);
  case CEIL:    printf ("ceil not implemented\n"); return new exprConst (0);
  case ABS:     return new exprAbs (nl2e (e -> L.e));
    //  case OPUMINUS:return new exprOpp (nl2e (e -> L.e -> L.e));
  case OPUMINUS:return new exprOpp (nl2e (e -> L.e));
  case OPIFnl:  printf ("ifnl not implemented\n"); return new exprConst (0);
  case OP_tanh: printf ("tanh not implemented\n"); return new exprConst (0);
  case OP_tan:  
    arg = nl2e (e -> L.e);
    return new exprDiv (new exprSin (arg), new exprCos (new exprClone (arg)));
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

  case OP_cos:     return new exprCos (nl2e (e -> L.e));
  case OP_atanh: printf ("atanh not implemented\n"); return new exprConst (0);
  case OP_atan2: printf ("atan2 not implemented\n"); return new exprConst (0);
  case OP_atan:  printf ("atan not implemented\n"); return new exprConst (0);
  case OP_asinh: printf ("asinh not implemented\n"); return new exprConst (0);
  case OP_asin:  printf ("asin not implemented\n"); return new exprConst (0);
  case OP_acosh: printf ("acosh not implemented\n"); return new exprConst (0);
  case OP_acos:  printf ("acos not implemented\n"); return new exprConst (0);
  case OPSUMLIST: {
    register int i=0;
    al = new expression * [(e->R.ep - e->L.ep)];
    for (ep = e->L.ep; ep < e->R.ep; ep++)
      al [i++] = nl2e (*ep);
    return new exprSum (al, i);
  }
  case OPintDIV: printf ("intdiv not implemented\n"); return new exprConst (0);
  case OPprecision: printf ("precision not implemented\n"); return new exprConst (0);
  case OPround:  printf ("round not implemented\n"); return new exprConst (0);
  case OPtrunc:  printf ("trunc not implemented\n"); return new exprConst (0);

  case OP1POW: return new exprPow (nl2e (e -> L.e), 
				   new exprConst (((expr_n *)e->R.e)->v));
  case OP2POW: return new exprPow (nl2e (e -> L.e), 
				   new exprConst (2));
  case OPCPOW: return new exprPow (new exprConst (((expr_n *)e->L.e)->v),
				   nl2e (e -> R.e));
  case OPFUNCALL: printf ("function call not implemented\n"); return new exprConst (0);
  case OPNUM:   return new exprConst (((expr_n *)e)->v);
  case OPPLTERM:  printf ("plterm not implemented\n"); return new exprConst (0);
  case OPIFSYM:   printf ("ifsym not implemented\n"); return new exprConst (0);
  case OPHOL:     printf ("hol not implemented\n"); return new exprConst (0);
  case OPVARVAL: return new exprClone (variables_ [e->a]);

  case -1:
  default: printf ("WARNING: unknown operator (address %x)\n", (long int) e -> op); 
    return new exprConst (0);
  }
}

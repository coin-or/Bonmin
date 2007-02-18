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


// warning for non-implemented functions -- return 0 constant expression

expression *notimpl (const std::string &fname) {
  std::cerr << "*** Warning: " << fname << " not implemented" << std::endl;
  return new exprConst (0);
}


// converts an AMPL expression (sub)tree into an expression* (sub)tree

expression *CouenneProblem::nl2e (expr2 *e) {

  expression **al;
  expr **ep;
  expression *arg;
  int j;

  switch (getOperator (e -> op)) {

  case OPPLUS:   return new exprSum (nl2e (e -> L.e), nl2e (e -> R.e));
  case OPMINUS:  return new exprSub (nl2e (e -> L.e), nl2e (e -> R.e -> L.e));
  case OPMULT:   return new exprMul (nl2e (e -> L.e), nl2e (e -> R.e));
  case OPDIV:    return new exprDiv (nl2e (e -> L.e), nl2e (e -> R.e));
  case OPREM:   return notimpl ("remainder");
  case OPPOW:    return new exprPow (nl2e (e -> L.e), nl2e (e -> R.e));
  case OPLESS:  return notimpl ("less");
  case MINLIST: return notimpl ("min");
  case MAXLIST: return notimpl ("max");
  case FLOOR:   return notimpl ("floor");
  case CEIL:    return notimpl ("ceil");
  case ABS:     return new exprAbs (nl2e (e -> L.e));
    //  case OPUMINUS:return new exprOpp (nl2e (e -> L.e -> L.e));
  case OPUMINUS:return new exprOpp (nl2e (e -> L.e));
  case OPIFnl:  return notimpl ("ifnl");
  case OP_tanh: return notimpl ("tanh");
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
  case OP_atanh: return notimpl ("atanh");
  case OP_atan2: return notimpl ("atan2");
  case OP_atan:  return notimpl ("atan");
  case OP_asinh: return notimpl ("asinh");
  case OP_asin:  return notimpl ("asin");
  case OP_acosh: return notimpl ("acosh");
  case OP_acos:  return notimpl ("acos");
  case OPSUMLIST: {
    register int i=0;
    al = new expression * [(e->R.ep - e->L.ep)];
    for (ep = e->L.ep; ep < e->R.ep; ep++)
      al [i++] = nl2e (*ep);
    return new exprSum (al, i);
  }
  case OPintDIV: return notimpl ("intdiv");
  case OPprecision: return notimpl ("precision");
  case OPround:  return notimpl ("round");
  case OPtrunc:  return notimpl ("trunc");

  case OP1POW: return new exprPow (nl2e (e -> L.e), 
				   new exprConst (((expr_n *)e->R.e)->v));
  case OP2POW: return new exprPow (nl2e (e -> L.e), 
				   new exprConst (2));
  case OPCPOW: return new exprPow (new exprConst (((expr_n *)e->L.e)->v),
				   nl2e (e -> R.e));
  case OPFUNCALL: return notimpl ("function call");
  case OPNUM:   return new exprConst (((expr_n *)e)->v);
  case OPPLTERM:  return notimpl ("plterm");
  case OPIFSYM:   return notimpl ("ifsym");
  case OPHOL:     return notimpl ("hol");
  case OPVARVAL:  
    // check if index is above number of variables
    // (should never happen, but it depends on some map [] in pfg_read.c)
    if ((j = ((expr2_v *) e) -> a) >= nVars ()) 
      for (int i=nVars(); i<=j; i++)
	addVariable (false);
    return new exprClone (variables_ [j]);

  case -1:
  default: printf ("WARNING: unknown operator (address %x)\n", (long int) e -> op); 
    return new exprConst (0);
  }
}

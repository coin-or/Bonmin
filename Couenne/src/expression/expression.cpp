/*
 * Name:    expression.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the expression class
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <string>
#include <string.h>

#include <CouenneCutGenerator.h>
#include <CouenneProblem.h>

#include <CouenneTypes.h>
#include <expression.h>
#include <exprAux.h>
#include <exprOp.h>
#include <exprUnary.h>
#include <exprVar.h>
#include <exprIVar.h>
#include <exprBound.h>

// static vectors for evaluation, see their description in
// expression.h

CouNumber  expression::stack [STACK_SIZE];
CouNumber *expression::sp = stack;

CouNumber *expression::variables_ = NULL;
CouNumber *expression::lbounds_   = NULL;
CouNumber *expression::ubounds_   = NULL;


// General N-ary function destructor

exprOp::~exprOp () {

  register expression *elem;

  if (arglist_) {
    for (register int i = nargs_; i--;)
      if ((elem = arglist_ [i]))
	delete elem;

    delete [] arglist_;
    arglist_ = NULL;
  }
}


// print expression

void exprOp::print (std::ostream      &out = std::cout, 
		    const std::string &op  = "unknown", 
		    enum pos           pos = PRE)        const 
{
  if (pos == PRE)
    out << op;

  out << "("; fflush (stdout);
  for (int i=0; i<nargs_; i++) {
    if (arglist_ [i])
      arglist_ [i] -> print (out); 
    fflush (stdout);
    if (i < nargs_ - 1) {
      if (pos == INSIDE) out << op;
      else               out << ",";
    }
    fflush (stdout);
  }
  out << ")";
  fflush (stdout);
}


// print unary expression

void exprUnary::print (std::ostream      &out = std::cout, 
		       const std::string &op = "unknown", 
		       enum pos           pos = PRE)       const 
{
  if (pos == PRE)  out << op;
  out << "(";
  argument_ -> print (out);
  out << ")";
  if (pos == POST) out << op;
}


// does this expression depend on variables in varlist?

bool exprOp::dependsOn (int *varlist = NULL, int n = 1) {

  for (register int i = nargs_; i--;)
    if (arglist_ [i] -> dependsOn (varlist, n))
      return true;

  return false;
}


// Get lower and upper bound of a generic expression

void expression::getBounds (expression *&lb, expression *&ub) {

  lb = new exprConst (- COUENNE_INFINITY);
  ub = new exprConst (  COUENNE_INFINITY);
}


// generate one cut for a constant

void exprConst::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			      OsiCuts &cs, const CouenneCutGenerator *cg) {

  OsiRowCut *cut;

  if ((cut = cg -> createCut (currValue_, 0, w -> Index (), 1.)))
    cs.insert (cut);
}


// name () -- a string value for each expression

#define MAX_NAME 10000

std::string Coutoa (CouNumber x) {
  char s [MAX_NAME];
  sprintf (s, "%f", x);
  return std::string (s);
}

std::string Indtoa (char prefix, int i) {
  char s [MAX_NAME];
  sprintf (s, "%c%d", prefix, i);
  return std::string (s);
}

const std::string exprVar::name   () const {return Indtoa ('x', varIndex_);}
const std::string exprIVar::name  () const {return Indtoa ('y', varIndex_);}
const std::string exprAux::name   () const {return Indtoa ('w', varIndex_);}

const std::string exprLowerBound::name () const {return Indtoa ('l', varIndex_);}
const std::string exprUpperBound::name () const {return Indtoa ('u', varIndex_);}

const std::string exprConst::name () const {return Coutoa (currValue_);}

const std::string exprCopy::name  () const {return copy_ -> Original () -> name ();}

const std::string exprUnary::name () const {return name ("?");}
const std::string exprOp::name    () const {return name ("?");}


/**
 *
 * The functions below could have been implemented much more easily as follows:
 *
 * const std::string exprUnary::name () const 
 *   {return "(" + argument_ -> name () + ")";}
 *
 * const std::string exprOp::   name () const {
 *  std::string args = "(" + arglist_ [0] -> name ();
 *  for (int i=1; i<nargs_; i++)
 *    args += "," + arglist_ [i] -> name ();
 *  return args + ")";
 * }
 *
 * But, because of a segfault occurred in some x86_64 machines, we have
 * decided to do it more plainly (and probably less efficiently).
 */

const std::string exprUnary::name (const std::string &op) const {

  register char *s = (char *) malloc (MAX_NAME * sizeof (char));
  sprintf (s, "%s(%s)", op.c_str (), argument_ -> name (). c_str ());

  s = (char *) realloc (s, (1 + strlen (s)) * sizeof (char));
  std::string ret (1 + strlen (s), ' ');
  for (int register i=strlen (s); i--;) // no need to set eos
    ret [i] = s [i]; 
  free (s);
  return ret;
}

const std::string exprOp::name (const std::string &op) const {

  register char *s = (char *) malloc (MAX_NAME * sizeof (char));
  sprintf (s, "%s(%s", op.c_str (), arglist_ [0] -> name (). c_str ());

  for (register int i=1; i<nargs_; i++) {
    strcat (s, ",");
    strcat (s, arglist_ [i] -> name (). c_str ());
  }

  strcat (s, ")");

  s = (char *) realloc (s, (1 + strlen (s)) * sizeof (char));
  std::string ret (1 + strlen (s), ' ');
  for (register int i=strlen (s); i--;) // no need to set eos
    ret [i] = s [i]; 
  free (s);
  return ret;
}

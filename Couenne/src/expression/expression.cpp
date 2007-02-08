/*
 * Name:    expression.C
 * Author:  Pietro Belotti
 * Purpose: methods of the expression class
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <sstream>

#include <CouenneTypes.h>
#include <expression.h>
#include <exprAux.h>
#include <exprOp.h>
#include <exprUnary.h>
#include <exprVar.h>
#include <exprIVar.h>
#include <exprBound.h>

#include <CouenneCutGenerator.h>
#include <CouenneProblem.h>

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

  out << "("; 
  for (int i=0; i<nargs_; i++) {
    if (arglist_ [i])
      arglist_ [i] -> print (out); 

    if (i < nargs_ - 1) {
      if (pos == INSIDE) out << op;
      else               out << ",";
    }
  }
  out << ")";
}


// print unary expression

void exprUnary::print (std::ostream      &out = std::cout, 
		       const std::string &op = "unknown", 
		       enum pos           pos = PRE)       const 
{
  if (pos == PRE)  out << op << " "; fflush (stdout);
  argument_ -> print (out);
  if (pos == POST) out << op; fflush (stdout);
}


// does this expression depend on variables in varlist?

bool exprOp::dependsOn (int *varlist = NULL, int n = 1) {

  for (register int i = nargs_; i--;)
    if (arglist_ [i] -> dependsOn (varlist, n))
      return true;

  return false;
}


// Get lower and upper bound of a generic expression

inline void expression::getBounds (expression *&lb, expression *&ub) {

  lb = new exprConst (- COUENNE_INFINITY);
  ub = new exprConst (  COUENNE_INFINITY);
}


// generate cuts for expression associated with this auxiliary

void exprAux::generateCuts (const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg) {

  //  printf ("Generate cut for "); 
  //  print (std::cout);  printf (" := ");
  //  image_ -> print (std::cout); printf("\n");

  int j = cs.sizeRowCuts ();
  image_ -> generateCuts (this, si, cs, cg);

  bool first=true;

  /*
  if (!(cg -> isFirst ()))
    for (;j < cs.sizeRowCuts ();j++) {

      if (first) {
	printf ("----------------Generated cut for "); 
	print (std::cout);  printf (" := ");
	image_ -> print (std::cout); printf("\n");
	first=false;
      }
      cs.rowCutPtr (j) -> print ();
    }
  */
}


// generate one cut for a constant

void exprConst::generateCuts (exprAux *w, const OsiSolverInterface &si, 
			      OsiCuts &cs, const CouenneCutGenerator *cg) {

  OsiRowCut *cut = new OsiRowCut;

  CouNumber *coeff = new CouNumber [1];  
  int       *index = new int       [1];  

  cut -> setLb (currValue_);
  cut -> setUb (currValue_);

  coeff [0] = 1; index [0] = w -> Index ();

  cs.insert (cut);
}


// name () -- a string value for each expression

std::string Coutoa (CouNumber x) {
  char s [50];
  sprintf (s, "%f", x);
  return std::string (s);
}

std::string Indtoa (int i) {
  char s [50];
  sprintf (s, "%d", i);
  return std::string (s);
}

const std::string exprVar::name   () const {return "x_" + Indtoa (varIndex_);}
const std::string exprIVar::name  () const {return "y_" + Indtoa (varIndex_);}
const std::string exprAux::name   () const {return "w_" + Indtoa (varIndex_);}

const std::string exprLowerBound::name () const {return "l_" + Indtoa (varIndex_);}
const std::string exprUpperBound::name () const {return "u_" + Indtoa (varIndex_);}

const std::string exprConst::name () const {return Coutoa (currValue_);}

const std::string exprCopy::name  () const {return copy_ -> Original () -> name ();}

const std::string exprUnary::name () const {return "(" + argument_ -> name() + ")";}

const std::string exprOp::name    () const {

  std::string args = "(" + arglist_ [0] -> name ();

  for (int i=1; i<nargs_; i++)
    args += "," + arglist_ [i] -> name ();

  return args + ")";
}

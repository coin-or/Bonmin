/*
 * Name:    problem.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the class CouenneProblem
 *
 * (C) Carnegie-Mellon University, 2006-07. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>

#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"

#include "CouenneTypes.hpp"

#include "expression.hpp"
#include "exprConst.hpp"
#include "exprGroup.hpp"
#include "exprClone.hpp"
#include "exprAux.hpp"

#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"

/// update value of variables, bounds

void CouenneProblem::update (const CouNumber *x, const CouNumber *l, const CouNumber *u, int n) {

  register int nvars = nVars ();

  // expand arrays if needed

  if (curnvars_ < nvars) {
    x_   = (CouNumber *) realloc (x_,  nvars * sizeof (CouNumber));
    lb_  = (CouNumber *) realloc (lb_, nvars * sizeof (CouNumber));
    ub_  = (CouNumber *) realloc (ub_, nvars * sizeof (CouNumber));

    curnvars_ = nvars;
  }

  // copy arrays (do not simply make x_ point to x)

  if (n < 0)
    n = nvars;

  {
    register int i;
    if (x) for (i=n; i--;) x_  [i] = x [i];
    if (l) for (i=n; i--;) lb_ [i] = l [i];
    if (u) for (i=n; i--;) ub_ [i] = u [i];
  }

  expression::update (x_, lb_, ub_);
}


/// initialize auxiliary variables from original variables in the
/// nonlinear problem

void CouenneProblem::initAuxs (CouNumber *x, 
			       CouNumber *l, 
			       CouNumber *u) {

  // update original variables only, that is, the first nVars ()
  // variables, as no auxiliaries exist yet
  update (x, l, u, nOrig_);

  // initially, auxiliary variables are unbounded, their bounds only
  // depending on their function

  for (register std::vector <exprVar *>::iterator i = variables_ . begin ();
       i != variables_ . end (); ++i)

    if (((*i) -> Type  () == AUX) &&    // this is an auxiliary
	((*i) -> Index () >= nOrig_)) { // and one that was not an original, originally...

      int index = (*i) -> Index ();
      lb_ [index] = - (ub_ [index] = COUENNE_INFINITY);
    }

  expression::update (x_, lb_, ub_);

  // only one loop is sufficient here, since auxiliary variable are
  // defined in such a way that w_i does NOT depend on w_j if i<j.

  Jnlst()->Printf(Ipopt::J_VECTOR, J_PROBLEM, "Initial bounds for aux (initAuxs):\n");

  for (int j=0, i=nVars (); i--; j++) {

    int ord = numbering_ [j];

    if (variables_ [ord] -> Type () == AUX) {

      exprAux *aux = dynamic_cast <exprAux *> (variables_ [ord]);

      Jnlst()->Printf(Ipopt::J_VECTOR, J_PROBLEM, "w_%04d [%10g,%10g] ", ord, lb_ [ord], ub_ [ord]);

      // set bounds 
      if ((lb_[ord] = CoinMax (lb_[ord], (*(aux -> Lb()))())) <= -COUENNE_INFINITY) lb_[ord]=-DBL_MAX;
      if ((ub_[ord] = CoinMin (ub_[ord], (*(aux -> Ub()))())) >=  COUENNE_INFINITY) ub_[ord]= DBL_MAX;
      //if ((lb_ [ord] = (*(aux -> Lb ())) ()) <= -COUENNE_INFINITY) lb_ [ord] = -DBL_MAX;
      //if ((ub_ [ord] = (*(aux -> Ub ())) ()) >=  COUENNE_INFINITY) ub_ [ord] =  DBL_MAX;

      Jnlst()->Printf(Ipopt::J_VECTOR, J_PROBLEM, " --> [%10g,%10g]\n", lb_ [ord], ub_ [ord]);

      x_ [ord] = CoinMax (lb_ [ord], CoinMin (ub_ [ord], (*(aux -> Image ())) ()));
    }
  }
}


/// get auxiliary variables from original variables in the nonlinear
/// problem

void CouenneProblem::getAuxs (CouNumber * x) {

  // save current addresses
  CouNumber 
    *xS = expression::Variables (),
    *lS = expression::Lbounds (),
    *uS = expression::Ubounds ();

  // temporarily make the expression arrays point to x (restore them
  // at the end of this function)
  expression::update (x, NULL, NULL);

  // set auxiliary w to f(x). This procedure is exact even though the
  // auxiliary variables have an incomplete image, i.e. they have been
  // decomposed previously, since they are updated with increasing
  // index.
  for (register int j = 0, i = nVars (); i--; j++) {

    exprVar *var = variables_ [numbering_ [j]];
    if (var -> Type () == AUX)
      x [var -> Index ()] = 
	//	(var -> Type () == AUX) ? 
	  (*(var -> Image ())) () 
	//	  (*var) ()
	;
  }

  // get the x and the bound vectors back to their previous state
  expression::update (xS, lS, uS);
}


/// fill obj vector with coefficient of the (linearized) obj function
/// (depends on sense of optimization -- invert if sense()==MAXIMIZE)

void CouenneProblem::fillObjCoeff (double *&obj) {

  // linearized objective can be an exprAux, an exprSub, an exprGroup,
  // or an exprSum. In the last two cases, the components are
  // variables or constants

  expression *body = objectives_ [0] -> Body ();
  int sense = objectives_ [0] -> Sense ();

  switch (body -> code ()) {

  case COU_EXPRVAR:   //
    obj [body -> Index ()] = (sense == MINIMIZE) ? 1 : -1;
    break;

  case COU_EXPRSUB: { // 

    expression **arglist = body -> ArgList ();

    obj [arglist [0] -> Index ()] = (sense == MINIMIZE) ?  1 : -1;
    obj [arglist [1] -> Index ()] = (sense == MINIMIZE) ? -1 :  1;

  } break;

  case COU_EXPRGROUP: { // 

    exprGroup *eg    = dynamic_cast <exprGroup *> (body);
    int       *index = eg -> getIndices ();
    CouNumber *coeff = eg -> getCoeffs  ();

    if (sense == MINIMIZE) while (*index >= 0) obj [*index++] =  *coeff++;
    else                   while (*index >= 0) obj [*index++] = -*coeff++;      
  } // no break, as exprGroup is derived from exprSum

  case COU_EXPRSUM: { // 

    expression **arglist = body -> ArgList ();

    for (int i = body -> nArgs (); i--;)
      switch ((arglist [i]) -> code ()) {

      case COU_EXPRCONST: 
	break;

      case COU_EXPRVAR: 
	obj [arglist [i] -> Index ()] = (sense == MINIMIZE) ? 1 : -1;
	break;

      case COU_EXPRMUL: {

	expression **mulArgList = arglist [i] -> ArgList ();
	int index = mulArgList [0] -> Index ();

	if (index >= 0) obj [index]                      = mulArgList [1] -> Value ();
	else            obj [mulArgList [1] -> Index ()] = mulArgList [0] -> Value ();
      }	break;

      default: 
	Jnlst()->Printf(Ipopt::J_ERROR, J_PROBLEM,
			"Couenne: invalid element of sum\nAborting\n");
	exit (-1);
      }
  } break;

  default:
    Jnlst()->Printf(Ipopt::J_WARNING, J_PROBLEM,
		    "### objective function not recognized\n");
    break;
  }
}


/// set cutoff from NLP solution
void CouenneProblem::setCutOff (CouNumber cutoff) 
{
  // AW: Should we use the value of the objecgive variable computed by 
  //     Couenne here?
  pcutoff_->setCutOff(cutoff + 1e-8 * fabs (cutoff));
  //  cutoff_ = cutoff + 1e-4 * fabs (cutoff);
} // tolerance needed to retain feasibility


/// Tell problem that auxiliary related to obj has a cutoff, to be
/// used in bound tightening
void CouenneProblem::installCutOff () {

  int indobj = objectives_ [0] -> Body () -> Index ();

  if (indobj >= 0) {

    //if (objectives_ [0] -> Sense () == MINIMIZE) 
    //{if (cutoff_ < ub_ [indobj]) ub_ [indobj] = cutoff_;}
    //else
    //{if (cutoff_ > lb_ [indobj]) lb_ [indobj] = cutoff_;}

    // all problem are assumed to be minimization
    double cutoff = pcutoff_->getCutOff();
    if (cutoff < ub_ [indobj])
      Jnlst()->Printf(Ipopt::J_DETAILED, J_PROBLEM,
		      "Installing cutoff %e for optimization variable index %d val = %e\n",
		      cutoff, indobj, ub_ [indobj]);
      ub_ [indobj] = cutoff;
  }
}

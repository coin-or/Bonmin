/*
 * Name:    problem.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the class CouenneProblem
 *
 * (C) Carnegie-Mellon University, 2006-08. 
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
#include "lqelems.hpp"
#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "lqelems.hpp"

const CouNumber SafeCutoff = COUENNE_EPS;

/// initialize auxiliary variables from original variables in the
/// nonlinear problem

void CouenneProblem::initAuxs () const {

  domain_.current () -> resize (nVars ());

  // initially, auxiliary variables are unbounded, their bounds only
  // depending on their function

  int nvars = nVars ();

  for (int i = 0; i < nvars; i++) {

    int index;

    if ((variables_ [i] -> Type  () == AUX) &&            // this is an auxiliary
	((index = variables_ [i] -> Index ()) >= nOrig_)) // and not an original, originally...
      //int index = variables_ [i] -> Index ();
      Lb (index) = - (Ub (index) = COUENNE_INFINITY);
  }

  // first initialize with values from constraints

  //Jnlst()->Printf(Ipopt::J_VECTOR, J_PROBLEM, "Initial bounds for aux (initAuxs):\n");

  for (std::vector <CouenneConstraint *>::const_iterator con = constraints_.begin ();
       con != constraints_.end (); ++con) {

    CouNumber
      lb = (*((*con) -> Lb ())) (),
      ub = (*((*con) -> Ub ())) ();

    int index = (*con) -> Body () -> Index ();

    assert (index >= 0);

    if ((Lb (index) = CoinMax (Lb (index), lb)) <= -COUENNE_INFINITY) Lb (index) = -COIN_DBL_MAX;
    if ((Ub (index) = CoinMin (Ub (index), ub)) >=  COUENNE_INFINITY) Ub (index) =  COIN_DBL_MAX;
  }

  // only one loop is sufficient here, since auxiliary variable are
  // defined in such a way that w_i does NOT depend on w_j if i<j.

  for (int j=0, i=nVars (); i--; j++) {

    int ord = numbering_ [j];

    if (variables_ [ord] -> Type () == AUX) {

      Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, "w_%04d [%10g,%10g] ", ord, Lb (ord), Ub (ord));

      CouNumber l, u;
      variables_ [ord] -> Image () -> getBounds (l, u);

      // set bounds 
      if ((Lb (ord) = CoinMax (Lb (ord), l)) <= -COUENNE_INFINITY) Lb (ord) = -COIN_DBL_MAX;
      if ((Ub (ord) = CoinMin (Ub (ord), u)) >=  COUENNE_INFINITY) Ub (ord) =  COIN_DBL_MAX;
      //if ((lb_ [ord] = (*(aux -> Lb ())) ()) <= -COUENNE_INFINITY) lb_ [ord] = -DBL_MAX;
      //if ((ub_ [ord] = (*(aux -> Ub ())) ()) >=  COUENNE_INFINITY) ub_ [ord] =  DBL_MAX;

      Jnlst () -> Printf (Ipopt::J_ALL, J_PROBLEM, " --> [%10g,%10g]\n", Lb (ord), Ub (ord));

      bool integer = variables_ [ord] -> isInteger ();

      if (integer) {
	Lb (ord) = ceil  (Lb (ord) - COUENNE_EPS);
	Ub (ord) = floor (Ub (ord) + COUENNE_EPS);
      }

      X (ord) = CoinMax (Lb (ord), CoinMin (Ub (ord), (*(variables_ [ord] -> Image ())) ()));
    }
  }
}


/// get auxiliary variables from original variables in the nonlinear
/// problem

void CouenneProblem::getAuxs (CouNumber * x) const {

  // set point at x, don't copy
  domain_.push (nVars (), x, domain_.lb (), domain_.ub (), false);

  // set auxiliary w to f(x). This procedure is exact even though the
  // auxiliary variables have an incomplete image, i.e. they have been
  // decomposed previously, since they are updated with increasing
  // index.
  for (int j=0, i=nVars (); i--; j++) {

    int index = numbering_ [j];
    exprVar *var = variables_ [index];

    CouNumber l, u;

    if (var -> Type () == AUX)
      var -> Image () -> getBounds (l,u);
    else {
      l = Lb (index);
      u = Ub (index);
    }

    if (var -> Type () == AUX)
      //x [index] = //not necessary, addresses of x and X are equal
      X (index) = 
	CoinMax (l, CoinMin (u, (*(var -> Image ())) ()));
  }

  domain_.pop ();
}


/// fill obj vector with coefficient of the (linearized) obj function
/// (depends on sense of optimization -- invert if sense()==MAXIMIZE)

void CouenneProblem::fillObjCoeff (double *&obj) {

  // linearized objective can be an exprAux, an exprSub, an exprGroup,
  // or an exprSum. In the last two cases, the components are
  // variables or constants

  expression *body = objectives_ [0] -> Body ();
  //int sense = objectives_ [0] -> Sense ();

  switch (body -> code ()) {

  case COU_EXPRVAR:   //
    obj [body -> Index ()] = 1; //(sense == MINIMIZE) ? 1 : -1;
    break;

  case COU_EXPRSUB: { // 

    expression **arglist = body -> ArgList ();

    obj [arglist [0] -> Index ()] =  1; //(sense == MINIMIZE) ?  1 : -1;
    obj [arglist [1] -> Index ()] = -1; //(sense == MINIMIZE) ? -1 :  1;

  } break;

  case COU_EXPRGROUP: { // 

    exprGroup *eg    = dynamic_cast <exprGroup *> (body);

    const exprGroup::lincoeff &lcoe = eg -> lcoeff ();

    //    if (sense == MINIMIZE) while (*index >= 0) obj [*index++] =  *coeff++;
    //    else                   while (*index >= 0) obj [*index++] = -*coeff++;      

    for (int n = lcoe.size (), i=0; n--; i++)
      //exprGroup::lincoeff::iterator el = lcoe.begin (); el != lcoe.end (); ++el)
      obj [lcoe [i]. first -> Index ()] = lcoe [i]. second;
    //(sense == MINIMIZE) ? 
    //(lcoe [i]. second) : 
    //-(lcoe [i]. second);

  } // no break, as exprGroup is derived from exprSum

  case COU_EXPRSUM: { // 

    expression **arglist = body -> ArgList ();

    for (int i = body -> nArgs (); i--;)
      switch ((arglist [i]) -> code ()) {

      case COU_EXPRCONST: 
	break;

      case COU_EXPRVAR: 
	obj [arglist [i] -> Index ()] = 1; //(sense == MINIMIZE) ? 1 : -1;
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

  case COU_EXPRCONST: break; // a constant objective

  default:
    Jnlst()->Printf(Ipopt::J_WARNING, J_PROBLEM,
		    "Couenne: warning, objective function not recognized\n");
    break;
  }
}


/// set cutoff from NLP solution
void CouenneProblem::setCutOff (CouNumber cutoff) const {

  int indobj = objectives_ [0] -> Body () -> Index ();

  // AW: Should we use the value of the objective variable computed by 
  //     Couenne here?
  if ((indobj >= 0) && (cutoff < pcutoff_ -> getCutOff () - COUENNE_EPS)) {

    Jnlst () -> Printf (Ipopt::J_DETAILED, J_PROBLEM,
			"Setting new cutoff %.10e for optimization variable index %d val = %.10e\n",
			cutoff, indobj,
			pcutoff_ -> getCutOff ());

    if (Var (indobj) -> isInteger ())
      pcutoff_    -> setCutOff (floor (cutoff + COUENNE_EPS));
    else pcutoff_ -> setCutOff (cutoff + SafeCutoff * fabs (1 + cutoff));
  }
} // tolerance needed to retain feasibility


/// Tell problem that auxiliary related to obj has a cutoff, to be
/// used in bound tightening
void CouenneProblem::installCutOff () {

  int indobj = objectives_ [0] -> Body () -> Index ();

  if (indobj >= 0) {

    // all problem are assumed to be minimization
    double cutoff = pcutoff_ -> getCutOff();

    if (cutoff < Ub (indobj))
      Ub (indobj) = cutoff;

  } else jnlst_ -> Printf (J_SUMMARY, J_PROBLEM, 
			   "Warning, could not install cutoff - negative objective index\n");
}


// clear all spurious variables pointers not referring to the variables_ vector
void CouenneProblem::realign () {

  // link variables to problem's domain
  for (std::vector <exprVar *>::iterator i = variables_.begin ();
       i != variables_.end (); ++i) {

    (*i) -> linkDomain (&domain_);
    (*i) -> realign (this);
    if ((*i) -> Type () == AUX)
      (*i) -> Image () -> realign (this);
  }

  // link variables to problem's domain
  for (std::vector <CouenneObjective *>::iterator i = objectives_.begin ();
       i != objectives_.end (); ++i) 
    (*i) -> Body () -> realign (this);


  // link variables to problem's domain
  for (std::vector <CouenneConstraint *>::iterator i = constraints_.begin ();
       i != constraints_.end (); ++i)
    (*i) -> Body () -> realign (this);
}

/// Add list of options to be read from file
void CouenneProblem::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> SetRegisteringCategory ("Couenne options", Bonmin::RegisteredOptions::CouenneCategory);

  roptions -> AddNumberOption
    ("art_cutoff",
     "Artificial cutoff",
     COIN_DBL_MAX,
     "Default value is infinity.");

  roptions -> AddNumberOption
    ("opt_window",
     "Window around known optimum",
     COIN_DBL_MAX,
     "Default value is infinity.");

  roptions -> AddNumberOption
    ("feas_tolerance",
     "Tolerance for constraints/auxiliary variables",
     feas_tolerance_default,
     "Default value is zero.");

  roptions -> AddStringOption2 
    ("feasibility_bt",
     "Feasibility-based (cheap) bound tightening",
     "yes",
     "no","",
     "yes","");

  roptions -> AddStringOption2 
    ("use_quadratic",
     "Use quadratic expressions and related exprQuad class",
     "no",
     "no","Use an auxiliary for each bilinear term",
     "yes","Create one only auxiliary for a quadrati expression");

  roptions -> AddStringOption2 
    ("optimality_bt",
     "Optimality-based (expensive) bound tightening",
     "yes",
     "no","",
     "yes","");

  roptions -> AddLowerBoundedIntegerOption
    ("log_num_obbt_per_level",
     "Specify the frequency (in terms of nodes) for optimality-based bound tightening.",
     -1,0,
     "\
If -1, apply at every node (expensive!). \
If 0, apply at root node only. \
If k>=0, apply with probability 2^(k - level), level being the current depth of the B&B tree.");

  roptions -> AddStringOption2 
    ("aggressive_fbbt",
     "Aggressive feasibility-based bound tightening (to use with NLP points)",
     "yes",
     "no","",
     "yes","");

  roptions -> AddLowerBoundedIntegerOption
    ("log_num_abt_per_level",
     "Specify the frequency (in terms of nodes) for aggressive bound tightening.",
     -1,2,
     "\
If -1, apply at every node (expensive!). \
If 0, apply at root node only. \
If k>=0, apply with probability 2^(k - level), level being the current depth of the B&B tree.");

  roptions -> AddNumberOption
    ("art_lower",
     "Artificial lower bound",
     -COIN_DBL_MAX,
     "Default value is -1.e50.");

  roptions -> AddStringOption3
    ("branching_object",
     "type of branching object for variable selection",
     "var_obj",
     "vt_obj",   "use Violation Transfer from Tawarmalani and Sahinidis",
     "var_obj",  "use one object for each variable",
     "expr_obj", "use one object for each nonlinear expression");
}

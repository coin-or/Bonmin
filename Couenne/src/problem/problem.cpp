/*
 * Name:    problem.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the class CouenneProblem
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>

//#include <sys/time.h>
#include <CoinHelperFunctions.hpp>
#include <CoinTime.hpp>

#include <CouenneTypes.h>
#include <expression.hpp>
#include <exprConst.hpp>
#include <exprClone.hpp>
#include <exprIVar.hpp>
#include <exprAux.hpp>
#include <exprMax.hpp>
#include <exprMin.hpp>

#include <CouenneProblem.hpp>
#include <CouenneProblemElem.hpp>


/// constructor
CouenneProblem::CouenneProblem (const struct ASL *asl):

  auxSet_    (NULL), 
  curnvars_  (-1),
  nIntVars_  (0),
  optimum_   (NULL),
  bestObj_   (COIN_DBL_MAX),
  quadIndex_ (NULL) {

  x_ = lb_ = ub_ = NULL; 

  if (!asl) return;

  double now = CoinCpuTime ();

  readnl (asl);

  if ((now = (CoinCpuTime () - now)) > 10.)
    printf ("Couenne: reading time %.3fs\n", now);

  now = CoinCpuTime ();
  //print (std::cout);
  //printf ("======================================\n");
  standardize ();

  fillQuadIndices ();

  if ((now = (CoinCpuTime () - now)) > 10.)
    printf ("Couenne: standardization time %.3fs\n", now);

  //readOptimum ("optimum-least.txt", optimum_, bestObj_, this);

  //print (std::cout);

  //problem_ -> writeMod ("extended-aw.mod", true);
  //problem_ -> writeMod ("extended-pb.mod", false);
}


/// clone problem

CouenneProblem *CouenneProblem::clone () const
  {return new CouenneProblem (*this);}


/// copy constructor

CouenneProblem::CouenneProblem (const CouenneProblem &p):
  x_        (NULL),
  lb_       (NULL),
  ub_       (NULL),
  curnvars_ (-1),
  nIntVars_ (p.nIntVars_),
  optimum_  (NULL),
  bestObj_  (p.bestObj_) {

  register int i;

  for (i=0; i < p.nObjs   (); i++) objectives_  . push_back (p.Obj   (i) -> clone ());
  for (i=0; i < p.nNLCons (); i++) constraints_ . push_back (p.NLCon (i) -> clone ());
  for (i=0; i < p.nVars   (); i++) variables_   . push_back (p.Var   (i) -> clone ());
  for (i=0; i < p.nAuxs   (); i++) auxiliaries_ . push_back (p.Aux   (i) -> clone ());

  if (p.optimum_) {
    optimum_ = (CouNumber *) malloc ((nVars () + nAuxs ()) * sizeof (CouNumber));

    for (i = nAuxs () + nVars (); i--;)
      optimum_ [i] = p.optimum_ [i];
  }

  update (p.X(), p.Lb(), p.Ub());
}


/// methods to add objective function

void CouenneProblem::addObjective (expression *newobj, const std::string &sense = "min") {
  objectives_ . push_back 
    (new Objective (newobj, (sense == "min") ? MINIMIZE : MAXIMIZE));
}


/// methods to add nonlinear constraints:

/// equality constraint
void CouenneProblem::addEQConstraint (expression *body, expression *rhs = NULL) {
  if (!rhs) rhs = new exprConst (0);
  constraints_ . push_back (new CouenneConstraint (body, rhs, new exprClone (rhs)));
}

/// "greater than" constraint
void CouenneProblem::addGEConstraint (expression *body, expression *rhs = NULL) {
  if (!rhs) rhs = new exprConst (0);
  constraints_ . push_back (new CouenneConstraint 
			    (body, rhs, new exprConst (1 + COUENNE_INFINITY)));
}

/// "smaller than" constraint
void CouenneProblem::addLEConstraint (expression *body, expression *rhs = NULL) {
  if (!rhs) rhs = new exprConst (0);
  constraints_ . push_back (new CouenneConstraint 
			    (body, new exprConst (- (1 + COUENNE_INFINITY)), rhs));
}

/// range constraint
void CouenneProblem::addRNGConstraint (expression *body, expression *lb=NULL, expression *ub=NULL) {
  if (!lb) lb = new exprConst (0);
  if (!ub) ub = new exprConst (0);
  constraints_ . push_back (new CouenneConstraint (body, lb, ub));
}


/// add variable to the problem -- check whether it is integer (isDiscrete)

expression *CouenneProblem::addVariable (bool isDiscrete) {

  exprVar *var = (isDiscrete) ? 
    (new exprIVar (variables_ . size () + auxiliaries_ . size ())) :
    (new exprVar  (variables_ . size () + auxiliaries_ . size ()));
  variables_ . push_back (var);

  if (isDiscrete) 
    nIntVars_++;

  return var;
}


/// add auxiliary variable and associate it with pointer to expression
/// given as argument

exprAux *CouenneProblem::addAuxiliary (expression *symbolic) {

  // check if image is already in the expression database auxSet_
  std::set <exprAux *, compExpr>::iterator i;

  // create new aux associated with that expression
  exprAux *w = new exprAux (symbolic, 
			    variables_ . size () + auxiliaries_ . size (), 
			    symbolic -> rank (this));

  // seek expression in the set
  if ((i = auxSet_ -> find (w)) == auxSet_ -> end ()) {

    // no such expression found in the set, create entry therein
    auxiliaries_ . push_back (w);
    auxSet_ -> insert (w);
  }
  else { // otherwise, just return the entry's pointer

    delete w;
    w = *i;
    w -> increaseMult ();
  }

  return w;
}


/// standardize all nonlinear objectives and constraints

void CouenneProblem::standardize () {

  // create expression set for binary search
  auxSet_ = new std::set <exprAux *, compExpr>;

  // standardize objectives
  for (std::vector <Objective *>::iterator i = objectives_.begin ();
       i != objectives_.end (); i++) {

    exprAux *aux = (*i) -> standardize (this);
    if (aux) (*i) -> Body (new exprClone (aux));
  }

  // standardize constraints
  for (std::vector <CouenneConstraint *>::iterator i = constraints_.begin ();
       i != constraints_.end (); i++) {

    exprAux *aux = (*i) -> standardize (this);
    if (aux) (*i) -> Body (new exprClone (aux));
  }

  delete auxSet_;

  int nTotVar = nVars() + nAuxs ();

  // reallocate space for enlarged set of variables
  x_  = (CouNumber *) realloc (x_,  nTotVar * sizeof (CouNumber));
  lb_ = (CouNumber *) realloc (lb_, nTotVar * sizeof (CouNumber));
  ub_ = (CouNumber *) realloc (ub_, nTotVar * sizeof (CouNumber));

  // make expression library point to new vectors
  expression::update (x_, lb_, ub_);

  for (int i=nVars (), j=0; j < nAuxs (); i++, j++) {

    // initial auxiliary bounds are infinite (they are later changed
    // through branching)

    lb_ [i] = -COUENNE_INFINITY;
    ub_ [i] =  COUENNE_INFINITY;

    // tighten them with propagated bounds
    auxiliaries_ [j] -> crossBounds ();

    // and evaluate them
    lb_ [i] = (*(auxiliaries_ [j] -> Lb    ())) ();
    ub_ [i] = (*(auxiliaries_ [j] -> Ub    ())) ();
    x_  [i] = (*(auxiliaries_ [j] -> Image ())) ();
  }
}


/// destroy problem components

CouenneProblem::~CouenneProblem () {

  if (x_) {
    free (x_);
    free (lb_);
    free (ub_);
  }

  // delete optimal solution (if any)
  if (optimum_)
    free (optimum_);

  // delete objectives
  for (std::vector <Objective *>::iterator i  = objectives_ . begin ();
       i != objectives_ . end (); i++)
    delete (*i);

  // delete constraints
  for (std::vector <CouenneConstraint *>::iterator i = constraints_ . begin ();
       i != constraints_ . end (); i++)
    delete (*i);

  // delete variables
  for (std::vector <exprVar *>::iterator i = variables_ . begin ();
       i != variables_ . end (); i++)
    delete (*i);

  // delete auxiliary variables
  for (std::vector <exprAux *>::iterator i = auxiliaries_ . begin ();
       i != auxiliaries_ . end (); i++)
    delete (*i);
}


/// update value of variables, bounds

void CouenneProblem::update (CouNumber *x, CouNumber *l, CouNumber *u, int n) {

  int nvars = nVars () + nAuxs ();

  // expand arrays if needed

  if (curnvars_ < nvars) {
    x_   = (CouNumber *) realloc (x_,  nvars * sizeof (CouNumber));
    lb_  = (CouNumber *) realloc (lb_, nvars * sizeof (CouNumber));
    ub_  = (CouNumber *) realloc (ub_, nvars * sizeof (CouNumber));

    curnvars_ = nvars;
  }

  // copy arrays (do not simply make x_ point to x)

  if (x) for (register int i = (n==-1) ? nvars : n; i--;)  x_  [i] = x [i];
  if (l) for (register int i = (n==-1) ? nvars : n; i--;)  lb_ [i] = l [i];
  if (u) for (register int i = (n==-1) ? nvars : n; i--;)  ub_ [i] = u [i];

  expression::update (x_, lb_, ub_);
}


/// initialize auxiliary variables from original variables in the
/// nonlinear problem

void CouenneProblem::initAuxs (CouNumber *x, 
			       CouNumber *l, 
			       CouNumber *u) {

  // update original variables only, that is, the first nVars ()
  // variables, as no auxiliaries exist yet
  update (x, l, u, nVars ());

  int nAux = nAuxs ();

  // initially, auxiliary variables are unbounded, their bounds only
  // depending on their function

  for (register int i=nVars (), j=nAux; j--; i++)
    lb_ [i] = - (ub_ [i] = COUENNE_INFINITY);

  expression::update (x_, lb_, ub_);

  // only one loop is sufficient here, since auxiliary variable are
  // defined in such a way that w_i does NOT depend on w_j if i<j.

  for (register int i = 0, j = nVars (); i < nAux; i++, j++) {

    exprAux *aux = Aux (i);

    x_ [j] = (*aux) ();

    // set bounds 
    if ((lb_ [j] = (*(aux -> Lb ())) ()) <= -COUENNE_INFINITY) lb_ [j] = -DBL_MAX;
    if ((ub_ [j] = (*(aux -> Ub ())) ()) >=  COUENNE_INFINITY) ub_ [j] =  DBL_MAX;
  }
}


/// get auxiliary variables from original variables in the nonlinear
/// problem

void CouenneProblem::getAuxs (CouNumber *x) {

  // temporarily make the expression arrays point to x (restore them
  // at the end of this function)
  expression::update (x, NULL, NULL);

  int nAux = nAuxs ();

  // set auxiliary w to f(x). This procedure is exact even though the
  // auxiliary variables have an incomplete image, i.e. they have been
  // decomposed previously, since they are updated with increasing
  // index.

  for (register int i = 0, j = nVars (); i < nAux; i++, j++)
    x [j] = (*(Aux (i) -> Image ())) ();

  // now get the x and the bound vectors back to their previous state
  expression::update (x_, lb_, ub_);
}

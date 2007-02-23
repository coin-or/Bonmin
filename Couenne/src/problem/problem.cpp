/*
 * Name:    problem.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the class CouenneProblem
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>

#include <CouenneTypes.h>
#include <expression.h>
#include <exprConst.h>
#include <exprClone.h>
#include <exprIVar.h>
#include <exprAux.h>
#include <exprMax.h>
#include <exprMin.h>

#include <CouenneProblem.h>
#include <CouenneProblemElem.h>

// clone problem

CouenneProblem *CouenneProblem::clone () const
  {return new CouenneProblem (*this);}


// copy constructor

CouenneProblem::CouenneProblem (const CouenneProblem &p) {

  register int i;

  for (i=0; i < p.nObjs   (); i++) objectives_  . push_back (p.Obj   (i) -> clone ());
  for (i=0; i < p.nNLCons (); i++) constraints_ . push_back (p.NLCon (i) -> clone ());
  for (i=0; i < p.nVars   (); i++) variables_   . push_back (p.Var   (i) -> clone ());
  for (i=0; i < p.nAuxs   (); i++) auxiliaries_ . push_back (p.Aux   (i) -> clone ());

  update (p.X(), p.Lb(), p.Ub());
}


// methods to add objective function

void CouenneProblem::addObjective (expression *newobj, const std::string &sense = "min") {
  objectives_ . push_back 
    (new Objective (newobj, (sense == "min") ? MINIMIZE : MAXIMIZE));
}


// methods to add nonlinear constraints:

// equality constraint
void CouenneProblem::addEQConstraint (expression *body, expression *rhs = NULL) {
  if (!rhs) rhs = new exprConst (0);
  constraints_ . push_back (new CouenneConstraint (body, rhs, new exprClone (rhs)));
}

// "greater than" constraint
void CouenneProblem::addGEConstraint (expression *body, expression *rhs = NULL) {
  if (!rhs) rhs = new exprConst (0);
  constraints_ . push_back (new CouenneConstraint 
			    (body, rhs, new exprConst (1 + COUENNE_INFINITY)));
}

// "smaller than" constraint
void CouenneProblem::addLEConstraint (expression *body, expression *rhs = NULL) {
  if (!rhs) rhs = new exprConst (0);
  constraints_ . push_back (new CouenneConstraint 
			    (body, new exprConst (- (1 + COUENNE_INFINITY)), rhs));
}

// range constraint
void CouenneProblem::addRNGConstraint (expression *body, expression *lb=NULL, expression *ub=NULL) {
  if (!lb) lb = new exprConst (0);
  if (!ub) ub = new exprConst (0);
  constraints_ . push_back (new CouenneConstraint (body, lb, ub));
}


// add variable to the problem

expression *CouenneProblem::addVariable (bool isDiscrete) {

  exprVar *var = (isDiscrete) ? 
    (new exprIVar (variables_ . size () + auxiliaries_ . size ())) :
    (new exprVar  (variables_ . size () + auxiliaries_ . size ()));
  variables_ . push_back (var);

  return var;
}


// add auxiliary variable and associate it with pointer to expression
// given as argument

exprAux *CouenneProblem::addAuxiliary (expression *symbolic) {

  // check if image is already in the expression database auxMap_

  exprAux *var;
  std::string key = symbolic -> name (); // string serving as a key to the map
  std::map <std::string, exprAux *>::iterator i;

  if ((i = auxMap_ -> find (key)) == auxMap_ -> end ()) {

    // no such expression has been found in the map, 

    // create entry in the map

    std::pair <std::string, exprAux *> newpair;
    newpair.first  = key;
    newpair.second = var = 
      // and corresponding auxiliary variable
      new exprAux (symbolic, variables_ . size () + auxiliaries_ . size ());
    auxiliaries_ . push_back (var);
    auxMap_ -> insert (newpair);
  }
  else var = (*i).second; // otherwise, just return the entry's
			  // auxiliary var. pointer

  return var;
}


// standardize all nonlinear objectives and constraints

void CouenneProblem::standardize () {

  // create expression map for binary search

  auxMap_ = new std::map <std::string, exprAux *>;

  // standardize objectives

  for (std::vector <Objective *>::iterator i = objectives_.begin ();
       i != objectives_.end (); i++) {

    exprAux *aux = (*i) -> standardize (this);

    if (aux)
      (*i) -> Body (new exprClone (aux));
  }

  // standardize constraints

  for (std::vector <CouenneConstraint *>::iterator i = constraints_.begin ();
       i != constraints_.end (); i++) {

    exprAux *aux = (*i) -> standardize (this);

    if (aux)
      (*i) -> Body (new exprClone (aux));
  }

  delete auxMap_;

  x_  = (CouNumber *) realloc (x_,  (nVars() + nAuxs ()) * sizeof (CouNumber));
  lb_ = (CouNumber *) realloc (lb_, (nVars() + nAuxs ()) * sizeof (CouNumber));
  ub_ = (CouNumber *) realloc (ub_, (nVars() + nAuxs ()) * sizeof (CouNumber));

  expression::update (x_, lb_, ub_);

  for (int i=0; i < nVars (); i++)
    (*(variables_ [i])) ();

  for (int i=nVars (), j=0; j < nAuxs (); i++, j++) {

    lb_ [i] = (*(auxiliaries_ [j] -> Lb    ())) ();
    ub_ [i] = (*(auxiliaries_ [j] -> Ub    ())) ();
    x_  [i] = (*(auxiliaries_ [j] -> Image ())) ();
  }
}


// destroy problem components

CouenneProblem::~CouenneProblem () {

  if (x_) {
    free (x_);
    free (lb_);
    free (ub_);
  }

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


// update value of variables, bounds

void CouenneProblem::update (CouNumber *x, CouNumber *l, CouNumber *u) {
  /*
  x_  = x;
  lb_ = l;
  ub_ = u;
  */

  static int curr_size = -1;

  int nvars = nVars () + nAuxs ();

  if (curr_size < nvars) {

    x_   = (CouNumber *) realloc (x_,  nvars * sizeof (CouNumber));
    lb_  = (CouNumber *) realloc (lb_, nvars * sizeof (CouNumber));
    ub_  = (CouNumber *) realloc (ub_, nvars * sizeof (CouNumber));

    curr_size = nvars;
  }

  for (register int i=nvars; i--;) {
    x_  [i] = x [i];
    lb_ [i] = l [i];
    ub_ [i] = u [i];
  }

  expression::update (x_, lb_, ub_);
}

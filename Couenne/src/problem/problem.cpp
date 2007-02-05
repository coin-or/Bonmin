/*
 * Name:    problem.C
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
  constraints_ . push_back (new CouenneConstraint (body, rhs, new exprCopy (rhs)));
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
  std::string key = symbolic -> name ();
  std::map <std::string, exprAux *>::iterator i;

  if ((i = auxMap_ -> find (key)) == auxMap_ -> end ()) {

    //    printf ("....... New exprAux!!! "); 
    //    symbolic -> print (std::cout);
    std::pair <std::string, exprAux *> newpair;
    newpair.first  = key;
    newpair.second = var = 
      new exprAux (symbolic, variables_ . size () + auxiliaries_ . size ());
    auxiliaries_ . push_back (var);
    auxMap_ -> insert (newpair);
    //    printf (" assigned to "); 
    //    var -> print (std::cout);
    //    printf ("\n");
  }
  else {
    //    symbolic -> print (std::cout);
    //    printf (" already assigned to "); 
    var = (*i).second;
    //    var -> print (std::cout);
    //    std::cout << "  (" << (*i).first << ")\n";
  }

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
      (*i) -> Body (aux);
  }

  // standardize constraints

  for (std::vector <CouenneConstraint *>::iterator i = constraints_.begin ();
       i != constraints_.end (); i++) {

    exprAux *aux = (*i) -> standardize (this);

    if (aux)
      (*i) -> Body (aux);
  }

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

  delete auxMap_;
}


// output content of the problem

void CouenneProblem::print (std::ostream &out = std::cout) {

  printf ("objectives:\n");
  for (std::vector <Objective *>::iterator i = objectives_.begin ();
       i != objectives_.end (); i++)
    (*i) -> print (out);

  printf ("constraints:\n");
  for (std::vector <CouenneConstraint *>::iterator i = constraints_.begin ();
       i != constraints_.end (); i++)
    (*i) -> print (out);

  printf ("auxiliaries:\n");
  for (std::vector <exprAux *>::iterator i = auxiliaries_.begin ();
       i != auxiliaries_.end (); i++) {
    (*i) -> print (out);             out << " := ";
    (*i) -> Image () -> print (out); out << std::endl;
  }
  printf ("end\n");
} 


// generate linear convexification of all nonlinear expressions
// associated with auxiliary variables
/*
void CouenneProblem::convexify () {

  expression ***coeff; // array of arrays of coefficients
  expression **rhs;    // array of arrays of right-hand sides
  int **indices;       // array of arrays of indices
  int *nterms;         // array with number of terms in each constraint
  enum con_sign *sign; // sign of all constraints added
  
  int n;               // number of tangents generated

  for (std::vector <exprAux *>::iterator aux = auxiliaries_.begin ();
       aux != auxiliaries_.end (); aux++) {

    // generate (set of) constraint to approximate function associated
    // with auxiliary variable (*i) from below, using linear ">="
    // constraints

    if ((n = (*aux) -> Image () -> lowerLinearHull (*aux, nterms, coeff, indices, rhs, sign))) {

      for (register int i=0; i<n; i++)
	addLinearConstraint (nterms [i], coeff [i], indices [i], rhs [i], sign [i]);

      delete [] coeff;     delete [] rhs;
      delete [] indices;   delete [] nterms;
      delete [] sign;
    }

    // same, but from above, using linear "<=" constraints

    if ((n = (*aux) -> Image () -> upperLinearHull (*aux, nterms, coeff, indices, rhs, sign))) {

      for (register int i=0; i<n; i++)
	addLinearConstraint (nterms [i], coeff [i], indices [i], rhs [i], sign [i]);

      delete [] coeff;     delete [] rhs;
      delete [] indices;   delete [] nterms;
      delete [] sign;
    }
  }

  // convexify those constraints that are originally linear 
  ********************** 
  for (std::vector <CouenneConstraint *>::iterator con = constraints_.begin ();
       con != constraints_.end (); con++)

    if ((*con) -> Body () -> Linearity () < QUADRATIC) {

      bool *var_dep = NULL;

      // this is a linear constraint: a sum, a difference, or a
      // multiple nesting of both, and as such has constant
      // derivatives. To avoid checking for dependency of all
      // variables, we build first a list of variables on which the
      // expression depends on.

      int n = (*con) -> Body () -> dependency (var_dep);

      expression **coeff   = new expression * [n];
      int         *indices = new int [n];

      int j=0;

      for (int i = nVar () + nAux (); i--;)

	if (var_dep [i]) {

	  coeff [j] = (*con) -> Body () -> differentiate (i);
	  coeff [j] -> simplify (); 
	  indices [j] = i;
	}

      addLinearConstraint (n, coeff, indices, 
			   (*con) -> Lb (), (*con) -> Ub (), (*con) -> sign ());

    }

      // generate (set of) constraint to approximate function associated
      // with auxiliary variable (*i) from below, using linear ">="
      // constraints

      if ((n = (*con) -> Body () -> lowerLinearHull (NULL, nterms, coeff, indices, rhs, sign))) {

	for (register int i=0; i<n; i++)
	  addLinearConstraint (nterms [i], coeff [i], indices [i], rhs [i], NULL, sign [i]);

	delete [] coeff;     delete [] rhs;
	delete [] indices;   delete [] nterms;
	delete [] sign;
      }

      // same, but from above, using linear "<=" constraints

      if ((n = (*con) -> Body () -> upperLinearHull (NULL, nterms, coeff, indices, rhs, sign))) {

	for (register int i=0; i<n; i++)
	  addLinearConstraint (nterms [i], coeff [i], indices [i], rhs [i], NULL, sign [i]);

	delete [] coeff;     delete [] rhs;
	delete [] indices;   delete [] nterms;
	delete [] sign;
      }
    }
  }
  ******************************************************************* /
}
*/


// destroy problem components

CouenneProblem::~CouenneProblem () {
  /*
  delete x_;
  delete lb_;
  delete ub_;
  */
  for (std::vector <CouenneConstraint *>::iterator i = constraints_ . begin ();
       i != constraints_ . end (); i++)
    delete (*i);

  for (std::vector <Objective *>::iterator i  = objectives_ . begin ();
       i != objectives_ . end (); i++)
    delete (*i);

  /*
  for (std::vector <LinearConstraint *>::iterator i = linearconstraints_ . begin ();
       i != linearconstraints_ . end (); i++)
    delete (*i);
  */
  /*
  for (std::vector <exprAux *>::iterator i = auxiliaries_ . begin ();
       i != auxiliaries_ . end (); i++)
    delete (*i);
  */
}


// update value of variables, bounds

void CouenneProblem::update (CouNumber *x, CouNumber *l, CouNumber *u) {

  int nvars = nVars () + nAuxs ();

  x_   = (CouNumber *) realloc (x_,  nvars * sizeof (CouNumber));
  lb_  = (CouNumber *) realloc (lb_, nvars * sizeof (CouNumber));
  ub_  = (CouNumber *) realloc (ub_, nvars * sizeof (CouNumber));

  register int i;

  for (i=nvars; i--;) {
    x_  [i] = x [i];
    lb_ [i] = l [i];
    ub_ [i] = u [i];
  }

  expression::update (x_, lb_, ub_);
}

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
#include <exprAux.h>
#include <exprMax.h>
#include <exprMin.h>

#include <CouenneProblem.h>
#include <CouenneProblemElem.h>

// clone problem

CouenneProblem *CouenneProblem::clone () const {

  return NULL;
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

expression *CouenneProblem::addVariable () {

  exprVar *var = new exprVar (variables_ . size () + auxiliaries_ . size ());
  variables_ . push_back (var);

  return var;
}


// add auxiliary variable and associate it with pointer to expression
// given as argument

exprAux *CouenneProblem::addAuxiliary (expression *added) {

  exprAux *var = new exprAux (added, variables_ . size () + auxiliaries_ . size ());
  auxiliaries_ . push_back (var);

  return var;
}


// standardize all nonlinear objectives and constraints

void CouenneProblem::standardize () {

  // standardize objectives

  for (std::vector <Objective *>::iterator i = objectives_.begin ();
       i != objectives_.end (); i++) {

      exprAux *aux = (*i) -> standardize (this);

      if (aux)
	(*i) -> Body (aux);
    }

  // same for constraints

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

  /*
  printf ("linear constraints:\n");
  for (std::vector <LinearConstraint *>::iterator i = linearconstraints_.begin ();
       i != linearconstraints_.end (); i++)
    (*i) -> print (out);
  */

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

// Allocate space in coeff, indices, rhs, and sign, for n constraint
// with number of coefficients given in nterms. Used in upper- and lowerLinearHull
/*
void allocateCon (int n, int *nterms,                     // input data
 		  expression ***& coeff, int **& indices, // allocated data
		  expression **& rhs, enum con_sign *& sign) {

  coeff   = new expression ** [n];
  indices = new int        *  [n];
  rhs     = new expression *  [n];
  sign    = new enum con_sign [n];

  while (n--) {
    register int nt = nterms [n];

    coeff   [n] = new expression * [nt];
    indices [n] = new int          [nt];
  }
}
*/

// destroy problem components

CouenneProblem::~CouenneProblem () {

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

  for (std::vector <exprAux *>::iterator i = auxiliaries_ . begin ();
       i != auxiliaries_ . end (); i++)
    delete (*i);
}

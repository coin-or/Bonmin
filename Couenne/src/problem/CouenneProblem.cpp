/*
 * Name:    CouenneProblem.cpp
 * Author:  Pietro Belotti
 * Purpose: main methods of the class CouenneProblem
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <vector>

#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"

#include "CouenneTypes.hpp"

#include "expression.hpp"
#include "exprConst.hpp"
#include "exprClone.hpp"
#include "exprIVar.hpp"
#include "exprAux.hpp"

#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "depGraph.hpp"

/// constructor
CouenneProblem::CouenneProblem (const struct ASL *asl,
				JnlstPtr jnlst):

  auxSet_    (NULL), 
  curnvars_  (-1),
  nIntVars_  (0),
  optimum_   (NULL),
  bestObj_   (COIN_DBL_MAX),
  quadIndex_ (NULL),
  commuted_  (NULL),
  numbering_ (NULL),
  ndefined_  (0),
  graph_     (NULL),
  nOrig_     (0),
  cutoff_    (COIN_DBL_MAX),
  jnlst_(jnlst) {

  x_ = lb_ = ub_ = NULL; 

  if (!asl) return;

  double now = CoinCpuTime ();

  readnl (asl);

  if ((now = (CoinCpuTime () - now)) > 10.)
    jnlst_->Printf (Ipopt::J_WARNING, J_PROBLEM,
		    "Couenne: reading time %.3fs\n", now);

  now = CoinCpuTime ();

  if (jnlst_->ProduceOutput(Ipopt::J_MOREVECTOR, J_PROBLEM)) {
    // we should put that also through the journalist
    print (std::cout);
    printf ("======================================\n");
  }


  // save -- for statistics purposes -- number of original
  // constraints. Some of them will be deleted as definition of
  // auxiliary variables.
  nOrigCons_ = constraints_. size ();

  // reformulation
  standardize ();

  // quadratic handling
  fillQuadIndices ();

  if ((now = (CoinCpuTime () - now)) > 10.)
    jnlst_->Printf(Ipopt::J_WARNING, J_CONVEXIFYING,
		   "Couenne: standardization time %.3fs\n", now);

  //  readOptimum ("nc.txt", optimum_, bestObj_, this);

  if (jnlst_->ProduceOutput(Ipopt::J_MOREVECTOR, J_PROBLEM)) {
    // We should route that also through the journalist
    print (std::cout);
  }

  //writeAMPL ("extended-aw.mod", true);
  //writeAMPL ("extended-pb.mod", false);
}


/// clone problem

CouenneProblem *CouenneProblem::clone () const
  {return new CouenneProblem (*this);}


/// copy constructor

CouenneProblem::CouenneProblem (const CouenneProblem &p):
  x_         (NULL),
  lb_        (NULL),
  ub_        (NULL),
  curnvars_  (-1),
  nIntVars_  (p.nIntVars_),
  optimum_   (NULL),
  bestObj_   (p.bestObj_),
  commuted_  (NULL),
  numbering_ (NULL),
  ndefined_  (p.ndefined_),
  graph_     (NULL),
  nOrig_     (p.nOrig_),
  nOrigCons_ (p.nOrigCons_),
  cutoff_    (p.cutoff_),
  jnlst_     (p.jnlst_) { // needed only in standardize (), unnecessary to update it

  // TODO: rebuild all lb_ and ub_ (needed for exprQuad)

  register int i;

  for (i=0; i < p.nObjs (); i++) objectives_  . push_back (p.Obj (i) -> clone ());
  for (i=0; i < p.nCons (); i++) constraints_ . push_back (p.Con (i) -> clone ());
  for (i=0; i < p.nVars (); i++) variables_   . push_back (p.Var (i) -> clone ());

  if (p.numbering_) {
    numbering_ = new int [i = nVars ()];
    while (i--)
      numbering_ [i] = p.numbering_ [i];
  }

  if (p.optimum_) {
    optimum_ = (CouNumber *) malloc (nVars () * sizeof (CouNumber));

    for (i = nVars (); i--;)
      optimum_ [i] = p.optimum_ [i];
  }

  update (p.X (), p.Lb (), p.Ub ());
}


/// destroy problem components

CouenneProblem::~CouenneProblem () {

  // free variables/bounds
  if (x_) {
    free (x_);
    free (lb_);
    free (ub_);
  }

  // delete optimal solution (if any)
  if (optimum_)
    free (optimum_);

  // delete objectives
  for (std::vector <CouenneObjective *>::iterator i  = objectives_ . begin ();
       i != objectives_ . end (); ++i)
    delete (*i);

  // delete constraints
  for (std::vector <CouenneConstraint *>::iterator i = constraints_ . begin ();
       i != constraints_ . end (); ++i)
    delete (*i);

  // delete variables
  for (std::vector <exprVar *>::iterator i = variables_ . begin ();
       i != variables_ . end (); ++i)
    delete (*i);

  // delete extra structures
  if (graph_)     delete    graph_;
  if (commuted_)  delete [] commuted_;
  if (numbering_) delete [] numbering_;
}


/// methods to add objective function

void CouenneProblem::addObjective (expression *newobj, const std::string &sense = "min") {
  objectives_ . push_back 
    (new CouenneObjective (newobj, (sense == "min") ? MINIMIZE : MAXIMIZE));
}


/// methods to add nonlinear constraints:

/// equality constraint
void CouenneProblem::addEQConstraint (expression *body, expression *rhs = NULL) {

  if (!rhs) rhs = new exprConst (0.);
  constraints_ . push_back (new CouenneConstraint (body, rhs, new exprClone (rhs)));
}

/// "greater than" constraint
void CouenneProblem::addGEConstraint (expression *body, expression *rhs = NULL) {
  if (!rhs) rhs = new exprConst (0.);
  constraints_ . push_back (new CouenneConstraint 
			    (body, rhs, new exprConst (COUENNE_INFINITY)));
}

/// "smaller than" constraint
void CouenneProblem::addLEConstraint (expression *body, expression *rhs = NULL) {
  if (!rhs) rhs = new exprConst (0.);
  constraints_ . push_back (new CouenneConstraint 
			    (body, new exprConst (- (COUENNE_INFINITY)), rhs));
}

/// range constraint
void CouenneProblem::addRNGConstraint (expression *body, expression *lb=NULL, expression *ub=NULL) {
  if (!lb) lb = new exprConst (0.);
  if (!ub) ub = new exprConst (0.);
  constraints_ . push_back (new CouenneConstraint (body, lb, ub));
}


/// add variable to the problem -- check whether it is integer (isDiscrete)

expression *CouenneProblem::addVariable (bool isDiscrete) {

  exprVar *var = (isDiscrete) ? 
    (new exprIVar (variables_ . size ())) :
    (new exprVar  (variables_ . size ()));

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
			    variables_ . size (), 
			    1 + symbolic -> rank (this));

  // seek expression in the set
  if ((i = auxSet_ -> find (w)) == auxSet_ -> end ()) {

    // no such expression found in the set, create entry therein
    variables_ . push_back (w);
    auxSet_ -> insert (w); // insert into repetition checking structure
    graph_  -> insert (w); // insert into acyclic structure

  } else {  // otherwise, just return the entry's pointer

    delete w;
    w = *i;
    (*i) -> increaseMult ();
  }

  return w;
}

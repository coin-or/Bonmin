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
#include "exprQuad.hpp"
#include "exprClone.hpp"
#include "exprIVar.hpp"
#include "exprAux.hpp"

#include "CouenneProblem.hpp"
#include "CouenneProblemElem.hpp"
#include "depGraph.hpp"
#include "lqelems.hpp"

/// constructor
CouenneProblem::CouenneProblem (const struct ASL *asl,
				Bonmin::BabSetupBase *base,
				JnlstPtr jnlst):
  problemName_ (""),
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
  pcutoff_   (new GlobalCutOff),
  created_pcutoff_ (true),
  doFBBT_    (false),
  doOBBT_    (false),
  doABT_     (false),
  logObbtLev_(0),
  jnlst_(jnlst),
  opt_window_ (1e50) {

  if (!asl) return;

  double now = CoinCpuTime ();

  // read problem from AMPL structure
  readnl (asl);

  if ((now = (CoinCpuTime () - now)) > 10.)
    jnlst_ -> Printf (Ipopt::J_WARNING, J_PROBLEM,
		      "Couenne: reading time %.3fs\n", now);

  now = CoinCpuTime ();

  // link initial variables to problem's domain
  for (std::vector <exprVar *>::iterator i = variables_.begin ();
       i != variables_.end (); ++i)
    (*i) -> linkDomain (&domain_);

  if (jnlst_->ProduceOutput(Ipopt::J_MOREVECTOR, J_PROBLEM))
    print (std::cout);

  // save -- for statistic purposes -- number of original
  // constraints. Some of them will be deleted as definition of
  // auxiliary variables.
  nOrigCons_ = constraints_. size ();

  jnlst_->Printf (Ipopt::J_SUMMARY, J_PROBLEM,
		  "Problem size before standarization: %d variables (%d integer) %d constraints.\n",
		  nOrig_, nIntVars (), nOrigCons_);

  // reformulation
  standardize ();

  // clear all spurious variables pointers not referring to the variables_ vector
  realign ();

  // quadratic handling
  fillQuadIndices ();

  if ((now = (CoinCpuTime () - now)) > 10.)
    jnlst_->Printf(Ipopt::J_WARNING, J_PROBLEM,
		   "Couenne: standardization time %.3fs\n", now);

  jnlst_->Printf(Ipopt::J_SUMMARY, J_PROBLEM,
		  "Problem size after standarization: %d variables (%d integer) %d constraints.\n",
		  nVars(), nIntVars(), nCons());

  if (jnlst_->ProduceOutput(Ipopt::J_MOREVECTOR, J_PROBLEM)) {
    // We should route that also through the journalist
    print (std::cout);
  }

  // give a value to all auxiliary variables
  initAuxs ();

  if (base) {

    std::string s;

    base -> options() -> GetStringValue ("feasibility_bt",     s, "couenne."); doFBBT_ = (s == "yes");
    base -> options() -> GetStringValue ("optimality_bt",      s, "couenne."); doOBBT_ = (s == "yes");
    base -> options() -> GetStringValue ("aggressive_fbbt",    s, "couenne."); doABT_  = (s == "yes");

    base -> options() -> GetIntegerValue ("log_num_obbt_per_level", logObbtLev_, "couenne.");
    base -> options() -> GetIntegerValue ("log_num_abt_per_level",  logAbtLev_,  "couenne.");

    CouNumber 
      art_cutoff =  2.e50,
      art_lower  = -2.e50;

    base -> options() -> GetNumericValue ("opt_window",  opt_window_, "couenne.");
    base -> options() -> GetNumericValue ("art_cutoff",  art_cutoff,  "couenne.");
    base -> options() -> GetNumericValue ("art_lower",   art_lower,   "couenne.");

    if (art_cutoff <  1.e50) setCutOff (art_cutoff);
    if (art_lower  > -1.e50) {
      int indobj = objectives_ [0] -> Body () -> Index ();
      if (indobj >= 0)
	domain_.lb (indobj) = art_lower;
    }
  }

  // check if optimal solution is available (for debug purposes)
  readOptimum ();

  //writeAMPL ("extended-aw.mod", true);
  //writeAMPL ("original.mod", false);
}


/// copy constructor

CouenneProblem::CouenneProblem (const CouenneProblem &p):
  problemName_ (p.problemName_),
  domain_    (p.domain_),
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
  pcutoff_   (p.pcutoff_),
  created_pcutoff_ (false),
  doFBBT_    (p. doFBBT_),
  doOBBT_    (p. doOBBT_),
  doABT_     (p. doABT_),
  logObbtLev_(p. logObbtLev_),
  logAbtLev_ (p. logAbtLev_),
  jnlst_     (p.jnlst_),
  opt_window_ (p.opt_window_) { // needed only in standardize (), unnecessary to update it

  for (int i=0; i < p.nVars (); i++)
    variables_ . push_back (NULL);

  //variables_ [0] = new exprVar (0, &domain_);

  for (int i=0; i < p.nVars (); i++) {
    int ind = p.numbering_ [i];
    //if (ind==0) delete variables_ [0];
    variables_ [ind] = p.Var (ind) -> clone (&domain_);

    //if (!(variables_ [ind])) printf ("var %d NULL\n", ind);
  }

  if (p.numbering_) {
    int i;
    numbering_ = new int [i = nVars ()];
    while (i--)
      numbering_ [i] = p.numbering_ [i];
  }

  for (int i=0; i < p.nObjs (); i++) objectives_  . push_back (p.Obj (i) -> clone (&domain_));
  for (int i=0; i < p.nCons (); i++) constraints_ . push_back (p.Con (i) -> clone (&domain_));

  if (p.optimum_) {
    optimum_ = (CouNumber *) malloc (nVars () * sizeof (CouNumber));

    for (int i = nVars (); i--;)
      optimum_ [i] = p.optimum_ [i];
  }

  // clear all spurious variables pointers not referring to the variables_ vector
  realign ();
}


/// Destructor

CouenneProblem::~CouenneProblem () {

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

  if (created_pcutoff_) delete pcutoff_;
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

expression *CouenneProblem::addVariable (bool isDiscrete, Domain *d) {

  exprVar *var = (isDiscrete) ? 
    (new exprIVar (variables_ . size (), d)) :
    (new exprVar  (variables_ . size (), d));

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
			    1 + symbolic -> rank (), 
			    exprAux::Unset, 
			    &domain_);
  //symbolic -> isInteger () ? exprAux::Integer : exprAux::Continuous);

  //  w -> linkDomain (&domain_);

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

/// Add list of options to be read from file
void CouenneProblem::registerOptions (Ipopt::SmartPtr <Bonmin::RegisteredOptions> roptions) {

  roptions -> SetRegisteringCategory ("Couenne options", Bonmin::RegisteredOptions::CouenneCategory);

  roptions -> AddNumberOption
    ("art_cutoff",
     "artificial cutoff",
     1.e50,
     "Default value is 1.e50.");

  roptions -> AddNumberOption
    ("opt_window",
     "window around known optimum",
     1.e5,
     "Default value is infinity.");

  roptions -> AddStringOption2 
    ("feasibility_bt",
     "Feasibility-based (cheap) bound tightening",
     "yes",
     "no","",
     "yes","");

  roptions -> AddStringOption2 
    ("optimality_bt",
     "optimality-based (expensive) bound tightening",
     "no",
     "no","",
     "yes","");

  roptions -> AddLowerBoundedIntegerOption
    ("log_num_obbt_per_level",
     "Specify the frequency (in terms of nodes) for optimality-based bound tightening.",
     -1,2,
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
     "artificial lower bound",
     -1.e50,
     "Default value is -1.e50.");
}


/// translates pair (indices, coefficients) into vector with pointers to variables
void CouenneProblem::indcoe2vector (int *indexL, 
				    CouNumber *coeff,
				    std::vector <std::pair <exprVar *, CouNumber> > &lcoeff) {
  for (int i=0; indexL [i] >= 0; i++)
    lcoeff.push_back (std::pair <exprVar *, CouNumber> (Var (indexL [i]), coeff [i]));
}


/// translates triplet (indicesI, indicesJ, coefficients) into vector with pointers to variables
void CouenneProblem::indcoe2vector (int *indexI,
				    int *indexJ,
				    CouNumber *coeff,
				    std::vector <quadElem> &qcoeff) {
  for (int i=0; indexI [i] >= 0; i++)
    qcoeff.push_back (quadElem (Var (indexI [i]), Var (indexJ [i]), coeff [i]));
}

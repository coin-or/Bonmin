/*
 * Name:    exprAux.cpp
 * Author:  Pietro Belotti
 * Purpose: methods of the class of auxiliary variables
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "exprAux.hpp"
#include "exprBound.hpp"
#include "exprMax.hpp"
#include "exprMin.hpp"
#include "CouenneCutGenerator.hpp"
#include "CouenneJournalist.hpp"

class CouenneCutGenerator;
class Domain;

//#define DEBUG

// auxiliary expression Constructor
exprAux::exprAux (expression *image, int index, int rank, enum intType isInteger, Domain *d): 

  exprVar       (index, d),
  image_        (image),
  rank_         (rank),
  multiplicity_ (1),
  integer_      (isInteger) {

  // do this later, in standardize()
  //  image_ -> getBounds (lb_, ub_);

  //  getBounds (lb_, ub_); // !!!!!!!!

  //  image_ -> getBounds (lb_, ub_);
  // getBounds (lb_, ub_);

  //  lb_ = new exprMax (new exprLowerBound (varIndex_), lb_);
  //  ub_ = new exprMin (new exprUpperBound (varIndex_), ub_);
  lb_ = new exprLowerBound (varIndex_, domain_);
  ub_ = new exprUpperBound (varIndex_, domain_);
}


/// Constructor to be used with standardize ([...], false)
exprAux::exprAux (expression *image, Domain *d):

  exprVar       (-1, d),
  image_        (image),
  lb_           (NULL),
  ub_           (NULL),
  rank_         (-1),
  multiplicity_ (0),
  integer_      (Unset) {}
//(image -> isInteger () ? Integer : Continuous)


/// Copy constructor
exprAux::exprAux (const exprAux &e, Domain *d):
  exprVar       (e.varIndex_, d), // variables ? (*variables) [0] -> domain () : NULL),
  image_        (e.image_ -> clone (d)),
  rank_         (e.rank_),
  multiplicity_ (e.multiplicity_),
  integer_      (e.integer_) {

  //  image_ -> getBounds (lb_, ub_);
  // getBounds (lb_, ub_);

  //  lb_ = new exprMax (new exprLowerBound (varIndex_), lb_);
  //  ub_ = new exprMin (new exprUpperBound (varIndex_), ub_);
  lb_ = new exprLowerBound (varIndex_, domain_);
  ub_ = new exprUpperBound (varIndex_, domain_);
}


/// Destructor
exprAux::~exprAux () {
  if (image_) 
    delete image_; 
  delete lb_;
  delete ub_;
}


/// Get lower and upper bound of an expression (if any)
void exprAux::getBounds (expression *&lb, expression *&ub) {

  // this replaces the previous 
  //
  //    image_ -> getBounds (lb0, ub0);
  //
  // which created large expression trees, now useless since all
  // auxiliaries are standardized.

  //  lb = lb_ -> clone ();//new exprLowerBound (varIndex_);
  //  ub = ub_ -> clone ();//new exprUpperBound (varIndex_);
  lb = new exprLowerBound (varIndex_, domain_);
  ub = new exprUpperBound (varIndex_, domain_);
}


/// set bounds depending on both branching rules and propagated
/// bounds. To be used after standardization
void exprAux::crossBounds () {

  expression *l0, *u0;

  image_ -> getBounds (l0, u0);

  //image_ -> getBounds (lb_, ub_);

  lb_ = new exprMax (lb_, l0);
  ub_ = new exprMin (ub_, u0);
}


/// I/O
void exprAux::print (std::ostream &out, bool descend) const {

  if (descend) 
    image_ -> print (out, descend);
  else {
    if (integer_) out << "z_"; // TODO: should be isInteger instead of
			       // integer_. Change all "isInteger()"
			       // to "isInteger() const"
    else          out << "w_";
    out << varIndex_;
  }
}


/// fill in the set with all indices of variables appearing in the
/// expression
int exprAux::DepList (std::set <int> &deplist, 
		      enum dig_type type) {

  if (type == ORIG_ONLY)   
    return image_ -> DepList (deplist, type);

  if (deplist.find (varIndex_) == deplist.end ())
    deplist.insert (varIndex_); 
  else return 0;

  if (type == STOP_AT_AUX) 
    return 1;

  return 1 + image_ -> DepList (deplist, type);
}


/// simplify
expression *exprAux::simplify () {

  if ((image_ -> Type () == AUX) || 
      (image_ -> Type () == VAR)) {

    --multiplicity_;
    expression *ret = image_;
    image_ = NULL;
    return ret;
  }

  return NULL;
}


// generate cuts for expression associated with this auxiliary

void exprAux::generateCuts (const OsiSolverInterface &si, 
			    OsiCuts &cs, const CouenneCutGenerator *cg, 
			    t_chg_bounds *chg, int,
			    CouNumber, CouNumber) {
  //#ifdef DEBUG
  static bool warned_large_coeff = false;
  int nrc = cs.sizeRowCuts (), ncc = cs.sizeColCuts ();
  //#endif

  /*
  if ((!(cg -> isFirst ())) && 
      ((l = domain_ -> lb (varIndex_)) > -COUENNE_INFINITY) &&
      ((u = domain_ -> ub (varIndex_)) <  COUENNE_INFINITY) &&
      (fabs (u-l) < COUENNE_EPS))
    cg -> createCut (cs, (l+u)/2., 0, varIndex_, 1.);
  else 
  */
  image_ -> generateCuts (this, si, cs, cg, chg);

  // check if cuts have coefficients, rhs too large or too small

  //#ifdef DEBUG

  if (cg -> Jnlst () -> ProduceOutput (Ipopt::J_MATRIX, J_CONVEXIFYING)) {
    if (warned_large_coeff)
      for (int jj=nrc; jj < cs.sizeRowCuts (); jj++) {

	OsiRowCut        *cut = cs.rowCutPtr (jj);
	CoinPackedVector  row = cut -> row ();

	int           n   = cut -> row (). getNumElements();
	const double *el  = row. getElements ();
	const int    *ind = row. getIndices ();
	double        rhs = cut -> rhs ();

	while (n--) {
	  if (fabs (el [n]) > COU_MAX_COEFF)  {
	    printf ("Couenne, warning: coefficient too large %g x%d: ", el [n], ind [n]);
	    cut -> print (); 
	    warned_large_coeff = true;
	    break;
	  }

	  if (fabs (rhs) > COU_MAX_COEFF) {
	    printf ("Couenne, warning: rhs too large (%g): ", rhs);
	    cut -> print ();
	    warned_large_coeff = true;
	    break;
	  }
	}
      }

    //  if (!(cg -> isFirst ())) 
    if (1 ||
	(nrc < cs.sizeRowCuts ()) || 
	(ncc < cs.sizeColCuts ()))
      {
	printf ("---------------- ConvCut:  "); 
	print (std::cout);  printf (" := ");
	image_ -> print (std::cout); 

	printf (" [%.7e,%.7e] <--- ", 
		domain_ -> lb (varIndex_), 
		domain_ -> ub (varIndex_));

	int index;
	if ((image_ -> Argument ()) && 
	    ((index = image_ -> Argument () -> Index ()) >= 0))
	  printf ("[%.7e,%.7e] ",
		  domain_ -> lb (index),
		  domain_ -> ub (index));
	else if (image_ -> ArgList ())
	  for (int i=0; i<image_ -> nArgs (); i++)
	    if ((index = image_ -> ArgList () [i] -> Index ()) >= 0)
	      printf ("[%.7e,%.7e] ", 
		      domain_ -> lb (index), 
		      domain_ -> ub (index));
	printf ("\n");

	for (int jj = nrc; jj < cs.sizeRowCuts (); jj++) cs.rowCutPtr (jj) -> print ();
	for (int jj = ncc; jj < cs.sizeColCuts (); jj++) cs.colCutPtr (jj) -> print ();
      }
  }
    //#endif

  //////////////////////////////////////////////////////////////

#if 0
  draw_cuts (cs, cg, nrc, this, image_);
#endif
}

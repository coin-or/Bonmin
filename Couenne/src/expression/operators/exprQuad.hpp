/*
 * Name:    exprQuad.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of quadratic expressions (constant+linear+quadratic)
 *
 * (C) Pietro Belotti 2007. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRQUAD_H
#define COUENNE_EXPRQUAD_H

#include <exprGroup.hpp>


///  class exprQuad, with constant, linear and quadratic terms

class exprQuad: public exprGroup {

 protected:

  /// Sparse implementation: given expression of the form sum_{i in N,
  /// j in N} a_{ij} x_i x_j, qindex0_ and qindex1_ contain
  /// respectively entries i and j for which a_{ij} is nonzero

  int       *qindexI_;
  int       *qindexJ_;
  CouNumber *qcoeff_;

 public:

  /// Constructor
  exprQuad  (CouNumber,                      /// constant part
	     int *, CouNumber *,             /// linear
	     int *, int *, CouNumber *,      /// quadratic
	     expression ** = NULL, int = 0); /// nonlinear

  /// copy constructor
  exprQuad (const exprQuad &src);

  /// destructor
  virtual ~exprQuad () {

    if (index_) {
      delete [] index_;
      delete [] coeff_;
    }

    if (qindexI_) {
      delete [] qindexI_;
      delete [] qindexJ_;
      delete [] qcoeff_;
    }
  }

  /// get indices and coefficients vectors of the quadratic part
  CouNumber *getQCoeffsI () {return qcoeff_;}
  int       *getQIndexI  () {return qindexI_;}
  int       *getQIndexJ  () {return qindexJ_;}

  /// cloning method
  virtual expression *clone () const
    {return new exprQuad (*this);}

  /// I/O
  virtual void print (std::ostream & = std::cout, bool = false, CouenneProblem * = NULL) const;

  /// function for the evaluation of the expression
  virtual CouNumber operator () ();

  /// differentiation
  virtual expression *differentiate (int index);

  /// simplification
  virtual expression *simplify ()
    {exprOp::simplify (); return NULL;}

  /// get a measure of "how linear" the expression is:
  virtual int Linearity ();/* {
    return ((qindexI_)                   ? QUADRATIC :
	    ((index_)                    ? LINEAR    :
	     ((fabs (c0_) > COUENNE_EPS) ? CONSTANT  : ZERO)));
	     }*/

  /// Get lower and upper bound of an expression (if any)
  virtual void getBounds (expression *&, expression *&) {}

  /// reduce expression in standard form, creating additional aux
  /// variables (and constraints)
  //  virtual exprAux *standardize (CouenneProblem *p);

  /// generate equality between *this and *w
  virtual void generateCuts (exprAux *w, const OsiSolverInterface &si, 
			     OsiCuts &cs, const CouenneCutGenerator *cg, 
			     t_chg_bounds * = NULL) {}

  /// only compare with people of the same kind
  virtual int compare (exprQuad &);

  /// code for comparisons
  virtual enum expr_type code () {return COU_EXPRQUAD;}

  /// used in rank-based branching variable choice
  virtual int rank (CouenneProblem *);

  /// return an index to the variable's argument that is better fixed
  /// in a branching rule for solving a nonconvexity gap
  virtual expression *getFixVar ();
};


/// compute sum of linear and nonlinear terms

inline CouNumber exprQuad::operator () () {

  register CouNumber  
     ret  = exprGroup::operator () (),
    *coe  = qcoeff_, 
    *vars = expression::Variables ();

  for (register int *qi = qindexI_, *qj = qindexJ_; *qi >= 0; )
    ret += *coe++ * vars [*qi++] * vars [*qj++];

  return (currValue_ = ret);
}

#endif

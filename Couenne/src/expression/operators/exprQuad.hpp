/**
 * Name:    exprQuad.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of quadratic expressions (= exprGroup +
 *          quadratic = constant + linear + [nonlinear] + quadratic)
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRQUAD_H
#define COUENNE_EXPRQUAD_H

#include <map>

#include "exprGroup.hpp"


/**  class exprQuad, with constant, linear and quadratic terms
 *
 *  It represents an expression of the form \f$a_0 + \sum_{i\in I} b_i
 *  x_i + x^T Q x + \sum_{i \in J} h_i (x)\f$, with \f$a_0 + \sum_{i\in
 *  I} b_i x_i\f$ an affine term, \f$x^T Q x\f$ a quadratic term, and
 *  a nonlinear sum \f$\sum_{i \in J} h_i (x)\f$. Standardization
 *  checks possible quadratic or linear terms in the latter and
 *  includes them in the former parts. 
 *
 *  If \f$h_i(x)\f$ is a product of two nonlinear, nonquadratic
 *  functions \f$h'(x)h''(x)\f$, two auxiliary variables
 *  \f$w'=f'(x)\f$ and \f$w''=h''(x)\f$ are created and the product
 *  \f$w'w''\f$ is included in the quadratic part of the exprQuad. If
 *  \f$h(x)\f$ nonquadratic, nonlinear function, an auxiliary variable
 *  \f$w=h(x)\f$ is created and included in the linear part.
 */

class exprQuad: public exprGroup {

 protected:

  /** \name Q matrix storage
    * Sparse implementation: given expression of the form \f$\sum_{i \in N,
    * j \in N} q_{ij} x_i x_j\f$, qindexI_ and qindexJ_ contain
    * respectively entries \f$i\f$ and \f$j\f$ for which \f$q_{ij}\f$ is nonzero
    * in \f$q_{ij} x_i x_j\f$: */
  /** @{ */

  int       *qindexI_; ///< the term \f$i\f$
  int       *qindexJ_; ///< the term \f$j\f$, (must be qindexJ_ [k] <= qindexI_ [k])
  CouNumber *qcoeff_;  ///< the term \f$a_{ij}\f$
  int        nqterms_; ///< number of non-zeroes in Q
  /** @} */

  /** \name Convexification data structures
   *
   *  These are filled by alphaConvexify through the method described
   *  in the LaGO paper by Nowak and Vigerske.
   */

  /// @{
  CouNumber *dCoeffLo_;  ///< diagonal coefficients of additional term for under convexfication
  CouNumber *dCoeffUp_;  ///< diagonal coefficients of additional term for over convexfication
  int       *dIndex_;    ///< indices of the \f$\alpha_i\f$'s in convexification
  int        nDiag_;     ///< number of elements in the above sparse vectors
  /// @}

 public:

  /** Constructor
   *
   * @param c0 constant part		     
   * @param li linear indices		     
   * @param lc linear coefficients		     
   * @param qi quadratic index \f$i\f$	     
   * @param qj quadratic index \f$j\f$	     
   * @param qc quadratic coefficient \f$a_{ij}\f$
   * @param al nonlinear arguments		     
   * @param n number of nonlinear arguments     
   */

  exprQuad  (CouNumber c0,
	     int *li,
	     CouNumber *lc,
	     int *qi,
	     int *qj,
	     CouNumber *qc,
	     expression **al = NULL,
	     int n = 0);

  /// Copy constructor
  exprQuad (const exprQuad &src);

  /// Destructor
  virtual ~exprQuad ();

  // get indices and coefficients vectors of the quadratic part
  CouNumber *getQCoeffs  () {return qcoeff_;}  ///< Get quadratic coefficients \f$q_{ij}\f$
  int       *getQIndexI  () {return qindexI_;} ///< Get index \f$i\f$ of the term \f$q_{ij}x_i x_j\f$
  int       *getQIndexJ  () {return qindexJ_;} ///< Get index \f$j\f$ of the term \f$q_{ij}x_i x_j\f$
  int        getnQTerms  () {return nqterms_;} ///< Get number of quadratic terms

  /// cloning method
  virtual expression *clone () const
    {return new exprQuad (*this);}

  /// Print expression to an iostream
  virtual void print (std::ostream & = std::cout, 
		      bool = false, CouenneProblem * = NULL) const;

  /// Function for the evaluation of the expression
  virtual CouNumber operator () ();

  /// Compute derivative of this expression with respect to variable
  /// whose index is passed as argument
  virtual expression *differentiate (int index);

  /// Simplify expression
  virtual expression *simplify ()
    {exprOp::simplify (); return NULL;}

  /// Get a measure of "how linear" the expression is
  virtual int Linearity () {
    int 
      lin  = exprSum::Linearity () >= NONLINEAR,
      lin2 = (qindexI_)                 ? QUADRATIC :
             (index_)                   ? LINEAR    :
	     (fabs (c0_) > COUENNE_EPS) ? CONSTANT  : ZERO;

    return ((lin > lin2) ? lin : lin2);
  }

  /// Get lower and upper bound of an expression (if any)
  virtual void getBounds (expression *&, expression *&);

  /// Generate cuts for the quadratic expression, which are supporting
  /// hyperplanes of the concave upper envelope and the convex lower
  /// envelope.
  virtual void generateCuts (exprAux *w, const OsiSolverInterface &si, 
			     OsiCuts &cs, const CouenneCutGenerator *cg, 
			     t_chg_bounds * = NULL, int = -1, 
			     CouNumber = -COUENNE_INFINITY, 
			     CouNumber =  COUENNE_INFINITY);

  /// Compute data for \f$\alpha\f$-convexification of a quadratic form
  /// (fills in dCoeff_ and dIndex_ for the convex underestimator)
  virtual void alphaConvexify (const OsiSolverInterface &);

  /** \function exprQuad::quadCuts 
   *
   * \brief Based on the information (dIndex_, dCoeffLo_, dCoeffUp_)
   * created/modified by alphaConvexify(), create convexification cuts
   * for this expression.
   *
   * The original constraint is :
   * \f[
   * \eta = a_0 + a^T x + x^T Q x
   * \f]
   * where \f$ \eta \f$ is the auxiliary corresponding to this
   * expression and \f$ w_j \f$ are the auxiliaries corresponding to
   * the other non-linear terms contained in the expression. 
   * 
   * The under-estimator of \f$ x^T Q x\f$ is given by \f[ x^T Q x +
   * \sum \lambda_{\min,i} (x_i - l_i ) (u_i - x_i ) \f] and its 
   * over-estimator is given by
   *
   * \f[ x^T Q x + \sum \lambda_{\max, i} (x_i - l_i ) (u_i - x_i )
   * \f] (where \f$ \lambda_{\min, i} = \frac{\lambda_{\min}}{w_i^2}
   * \f$, \f$ \lambda_{\max, i} = \frac{\lambda_{\max}}{w_i^2} \f$,
   * and \f$w_i = u_i - l_i\f$).
   *
   * Let \f$ \tilde a_0(\lambda)\f$, \f$ \tilde a(\lambda) \f$ and
   * \f$ \tilde Q (\lambda) \f$ be
   *
   * \f[ \tilde a_0(\lambda) = a_0 - \sum_{i = 1}^n \lambda_i l_i u_i \f]
   *
   * \f[ \tilde a(\lambda) = a + \left[ \begin{array}{c} \lambda_1
   * (u_1 + l_1) \\ \vdots \\ \lambda_n (u_n + l_n) \end{array}
   * \right], \f]
   *
   * \f[ \tilde Q(\lambda) = Q - \left( \begin{array}{ccc}
   * {\lambda_1} & & 0 \\ & \ddots & \\ 0 & & \lambda_n \end{array}
   * \right). \f]
   *
   * The convex relaxation of the initial constraint is then given by
   * the two constraints
   *
   * \f[ \eta \geq \tilde a_0(\lambda_{\min}) + \tilde
   * a(\lambda_{\min})^T x + x^T \tilde Q(\lambda_{\min}) x \f]
   *
   * \f[ \eta \leq \tilde a_0(\lambda_{\max}) + \tilde
   * a(\lambda_{\max})^T x + x^T \tilde Q(\lambda_{\max}) x \f]
   *  
   * The cut is computed as follow. Let \f$ (x^*, \eta^*) \f$ be
   * the solution at hand. The two outer-approximation cuts are:
   *
   * \f[ \eta \geq \tilde a_0(\lambda_{\min}) + \tilde
   * a(\lambda_{\min})^T x + {x^*}^T \tilde Q(\lambda_{\min}) (2x -
   * x^*) \f]
   *
   * and
   *
   * \f[ \eta \leq \tilde a_0(\lambda_{\max}) + \tilde
   * a(\lambda_{\max})^T x + {x^*}^T \tilde Q(\lambda_{\max}) (2x -
   * x^*) \f]
   *
   * grouping coefficients, we get:
   *
   * \f[ {x^*}^T \tilde Q(\lambda_{\min}) x^* - \tilde
   * a_0(\lambda_{\min}) \geq ( a(\lambda_{\min}) + 2
   * \tilde Q(\lambda_{\min} ) x^*)^T x - \eta \f]
   *
   * and
   *
   * \f[ {x^*}^T \tilde Q(\lambda_{\max}) x^* - \tilde
   * a_0(\lambda_{\max}) \leq ( a(\lambda_{\max}) + 2
   * \tilde Q (\lambda_{\max}) x^* )^T x - \eta \f]
   */

  void quadCuts (exprAux *w, OsiCuts & cs, const CouenneCutGenerator * cg);

  /// Compare two exprQuad
  virtual int compare (exprQuad &);

  /// Code for comparisons
  virtual enum expr_type code () {return COU_EXPRQUAD;}

  /// Used in rank-based branching variable choice
  virtual int rank (CouenneProblem *);

  /// is this expression integer? TODO: see exprGroup...
  virtual inline bool isInteger ()
    {return false;}

  /// fill in the set with all indices of variables appearing in the
  /// expression
  virtual int DepList (std::set <int> &deplist, 
		       enum dig_type type = ORIG_ONLY,
		       CouenneProblem *p = NULL);

  /// Return an index to the variable's argument that is better fixed
  /// in a branching rule for solving a nonconvexity gap. For this
  /// expression, always return NULL as selectBranch() will return the
  /// correct variable
  virtual expression *getFixVar ();

  /// Set up branching object by evaluating many branching points for
  /// each expression's arguments
  virtual CouNumber selectBranch (const CouenneObject *, const OsiBranchingInformation *,
				  int &, double * &, int &);

  /// Fill dependence set of the expression associated with this
  /// auxiliary variable
  virtual void fillDepSet (std::set <DepNode *, compNode> *dep, DepGraph *g);

  /// implied bound processing
  virtual bool impliedBound (int, CouNumber *, CouNumber *, t_chg_bounds *);

  /// create dIndex_ based on occurrences in qindexI_ and qindexJ_
  void make_dIndex (int, int *);
};


/// Compute sum of linear and nonlinear terms

inline CouNumber exprQuad::operator () () {

  register CouNumber
     ret  = exprGroup::operator () (), // compute non-quadratic part (linear+constant)
    *coe  = qcoeff_, 
    *vars = expression::Variables ();

  for (register int *qi = qindexI_, 
	            *qj = qindexJ_, 
	 i = nqterms_; i--; qi++, qj++) {

    register double term = 
      *coe++ * vars [*qi] * vars [*qj];

    ret += (*qi == *qj) ? term : 2 * term;
  }

  return ret;
}


/// Insert a pair <int,CouNumber> into a map for linear terms
void linsert (std::map <int, CouNumber> &lmap, 
	      int index, CouNumber coe);

/// Insert a pair <<int,int>,CouNumber> into a map for quadratic terms
void qinsert (std::map <std::pair <int, int>, CouNumber> &map, 
	      int indI, int indJ, CouNumber coe);

#endif

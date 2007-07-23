/*
 * Name:    exprQuad.hpp
 * Author:  Pietro Belotti
 * Purpose: definition of quadratic expressions (= exprGroup +
 *          quadratic = constant + linear + [nonlinear] + quadratic)
 *
 * (C) Pietro Belotti 2007. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRQUAD_H
#define COUENNE_EXPRQUAD_H

#include <exprGroup.hpp>


/**  class exprQuad, with constant, linear and quadratic terms
 *
 *
 *
 *
 *
 *
 *
 *
 *
 */

class exprQuad: public exprGroup {

 protected:

  /** \name Q matrix storage
    * Sparse implementation: given expression of the form sum_{i in N,
    * j in N} a_{ij} x_i x_j, qindexI_ and qindexJ_ contain
    * respectively entries i and j for which a_{ij} is nonzero
    * in a_ij x_i x_j: */
  /** @{ */

  int       *qindexI_; ///< the term i
  int       *qindexJ_; ///< the term j, (must be qindexJ_ [k] <= qindexI_ [k])
  CouNumber *qcoeff_;  ///< the term a_ij
  int        nqterms_; ///< number of non-zeroes in Q
  /** @} */

  CouNumber *dCoeffLo_;  ///< diagonal coefficients of additional term for under convexfication
  CouNumber *dCoeffUp_;  ///< diagonal coefficients of additional term for over convexfication
  int       *dIndex_;    ///< and indices (this is a sparse vector)
  int        nDiag_;     ///< number of elements in the above sparse vectors

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
  virtual ~exprQuad () {

    if (qindexI_) {
      delete [] qindexI_;
      delete [] qindexJ_;
      delete [] qcoeff_;
    }

    if (dIndex_) {
      delete [] dIndex_;
      delete [] dCoeffLo_;
      delete [] dCoeffUp_;
    }
  }

  /// get indices and coefficients vectors of the quadratic part
  CouNumber *getQCoeffs  () {return qcoeff_;}
  int       *getQIndexI  () {return qindexI_;}
  int       *getQIndexJ  () {return qindexJ_;}
  int        getnQTerms  () {return nqterms_;}

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

  /// generate equality between *this and *w
  virtual void generateCuts (exprAux *w, const OsiSolverInterface &si, 
			     OsiCuts &cs, const CouenneCutGenerator *cg, 
			     t_chg_bounds * = NULL, int = -1, 
			     CouNumber = -COUENNE_INFINITY, 
			     CouNumber =  COUENNE_INFINITY);

  /// [Stefan] fills in dCoeff_ and dIndex_ for the convex
  /// underestimator of this expression
  virtual void alphaConvexify (const OsiSolverInterface &);

  /** \function exprQuad::quadCuts 
    * \brief Based on the information (dIndex_, dCoeffLo_, dCoeffUp_)
    * created/modified by alphaConvexify(), create convexification cuts
    * for this expression.
    *
    * The original constraint is :
    * \f[
    * \eta = a_0 + a^T x + x^T Q x + \sum w_j
    * \f]
    * where \f$ \eta \f$ is the auxiliary corresponding to this
    * expression and \f$ w_j \f$ are the auxiliaries corresponding to
    * the other non-linear terms contained in the expression. (I don't
    * think that it is assumed anywhere in the function that Q only
    * involves original variables of the problem).
    * 
    * The under-estimator of \f$ x^T Q x\f$ is given by \f$ x^T Q x +
    * \sum \lambda_{\min,i} (x_i - l_i ) (u_i - x_i )} \f$ and its
    * over-estimator is given by \f$ Q - \sum \lambda_{\max, i} (x_i -
    * l_i ) (u_i - x_i ) \f$ (where \f$ \lambda_{\min, i} =
    * \frac{\lambda_{\min}}{w_i^2} \f$ and \f$ \lambda_{\max, i} =
    * \frac{\lambda_{\max}}{w_i^2} \f$).
    *
    * Let \f$ \tilde a_0(\lambda)\f$, \f$ \tilde a(\lambda) \f$ and
    * \f$ \tilde Q_(\lambda) \f$ be
    *
    * \f[ \tilde a_0(\lambda) = a_0 - \sum_{i = 1}^n \lambda_i l_i u_i \f]
    *
    * \f[ \tilde a(\lambda) = a + \left[ \begin{array}{c} \lambda_1
    * (u_1 + l_1) \\ \vdots \\ \lambda_n (u_n + l_n) \end{array}
    * \right], \f]
    *
    * \f[ \tilde Q(\lambda) = Q - \left( \begin{array}{ccc}
    * {\lambda_1} & 0 \\ & \ddots & \\ 0 & & \lambda_n \end{array}
    * \right). \f]
    *
    * The convex relaxation of the initial constraint is then given by
    * the two constraints
    *
    * \f[ \eta \geq \tilde a_0(\lambda_{\min}) + \tilde
    * a(\lambda_{\min})^T x + x^T \tilde Q(\lambda_{\min}) x + \sum
    * z_j \f]
    *
    * \f[ \eta \leq \tilde a_0(- \lambda_{\max}) + \tilde
    * a(-\lambda_{\max})^T x + x^T \tilde Q(-\lambda_{\max}) x + \sum
    * z_j \f]
    *  
    * The cut is computed as follow. Let \f$ (x^*, z^*, \eta^*) \f$ be
    * the solution at hand. The two outer-approximation cuts are:
    *
    * \f[ \eta \geq \tilde a_0(\lambda_{\min}) + \tilde
    * a(\lambda_{\min})^T x + {x^*}^T \tilde Q(\lambda_{\min}) (2x -
    * x^*) + \sum z_j \f]
    *
    * and
    *
    * \f[ \eta \leq \tilde a_0(-\lambda_{\max}) + \tilde
    * a(-\lambda_{\max})^T x + {x^*}^T \tilde Q(-\lambda_{\max}) (2x -
    * x^*) + \sum z_j \f]
    *
    * grouping coefficients, we get:
    *
    * \f[ {x^*}^T \tilde Q(\lambda_{\min}) x^* - \tilde
    * a_0(\lambda_{\min}) \geq ( a(\lambda_{\min}) + 2
    * Q(\lambda_{\min} ) x^*)^T x + \sum z_j - \eta \f]
    *
    * and
    *
    * \f[ {x^*}^T \tilde Q(-\lambda_{\max}) x^* - \tilde
    * a_0(-\lambda_{\max}) \leq ( a(-\lambda_{\max}) + 2
    * Q(-\lambda_{\max}) x^* )^T x + \sum z_j - \eta \f]
    */

  void quadCuts (exprAux *w, OsiCuts & cs, const CouenneCutGenerator * cg);

  /// only compare with people of the same kind
  virtual int compare (exprQuad &);

  /// code for comparisons
  virtual enum expr_type code () {return COU_EXPRQUAD;}

  /// used in rank-based branching variable choice
  virtual int rank (CouenneProblem *);

  /// return an index to the variable's argument that is better fixed
  /// in a branching rule for solving a nonconvexity gap
  virtual expression *getFixVar ();

  /// set up branching object by evaluating many branching points for
  /// each expression's arguments
  virtual CouNumber selectBranch (expression *, const OsiBranchingInformation *,
				  int &, double * &, int &);

  /// fill dependence set of the expression associated with this
  /// auxiliary variable
  virtual void fillDepSet (std::set <DepNode *, compNode> *dep, DepGraph *g);
};


/// compute sum of linear and nonlinear terms

inline CouNumber exprQuad::operator () () {

  register CouNumber  
     ret  = exprGroup::operator () (),
    *coe  = qcoeff_, 
    *vars = expression::Variables ();

  for (register int *qi = qindexI_, *qj = qindexJ_, i = nqterms_; i--;)

    ret += (*qi == *qj) ?
      (*coe++ *     vars [*qi++] * vars [*qj++]) :
      (*coe++ * 2 * vars [*qi++] * vars [*qj++]);

  return (currValue_ = ret);
}

#endif

// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Author:   Andreas Waechter                 IBM    2007-03-29

#ifndef __BONCUTSTRENGTHENER_HPP__
#define __BONCUTSTRENGTHENER_HPP__

#include "BonTMINLP.hpp"
#include "CoinPackedVector.hpp"
#include "BonTNLPSolver.hpp"

namespace Bonmin
{
  using namespace Ipopt;

  enum CutStrengtheningType{
    CS_None=0,
    CS_StrengthenedGlobal=1,
    CS_UnstrengthenedGlobal_StrengthenedLocal=2,
    CS_StrengthenedGlobal_StrengthenedLocal=3
  };

  enum DisjunctiveCutType{
    DC_None=0,
    DC_MostFractional=1
  };

  /** Class for strengthening OA cuts, and generating additional ones.
   */
  class CutStrengthener: public ReferencedObject
  {
    /** Class implementing the TNLP for strengthening one cut.  We
     *  assume that the cut has a lower bound. */
    class StrengtheningTNLP: public TNLP {
    public:
      /** Contructor */
      StrengtheningTNLP(SmartPtr<TMINLP> tminlp,
			const CoinPackedVector& cut,
			bool lower_bound,
			Index n,
			const Number* starting_point,
			const double* x_l_orig,
			const double* x_u_orig,
			Index constr_index,
			Index nvar_constr /** Number of variables in constraint */,
			const Index* jCol);

      /** Destructor */
      ~StrengtheningTNLP();

      /**@name Overloaded from TNLP */
      //@{
      /** Method to return some info about the nlp */
      virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
				Index& nnz_h_lag, IndexStyleEnum& index_style);

      /** Method to return the bounds for my problem */
      virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
				   Index m, Number* g_l, Number* g_u);

      /** Method to return the starting point for the algorithm */
      virtual bool get_starting_point(Index n, bool init_x, Number* x,
				      bool init_z, Number* z_L, Number* z_U,
				      Index m, bool init_lambda,
				      Number* lambda);

      /** Method to return the objective value */
      virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

      /** Method to return the gradient of the objective */
      virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

      /** Method to return the constraint residuals */
      virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

      /** Method to return:
       *   1) The structure of the jacobian (if "values" is NULL)
       *   2) The values of the jacobian (if "values" is not NULL)
       */
      virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
			      Index m, Index nele_jac, Index* iRow, Index *jCol,
			      Number* values);

      /** Method to return:
       *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
       *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
       */
      virtual bool eval_h(Index n, const Number* x, bool new_x,
			  Number obj_factor, Index m, const Number* lambda,
			  bool new_lambda, Index nele_hess, Index* iRow,
			  Index* jCol, Number* values);

      //@}
      /** @name Solution Methods */
      //@{
      /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
      virtual void finalize_solution(SolverReturn status,
				     Index n, const Number* x, const Number* z_L, const Number* z_U,
				     Index m, const Number* g, const Number* lambda,
				     Number obj_value,
				     const IpoptData* ip_data,
				     IpoptCalculatedQuantities* ip_cq);
      //@}

      /** Method for asking for the strengthened bound. */
      Number StrengthenedBound() const;
      
    private:
      /**@name Methods to block default compiler methods. */
      //@{
      StrengtheningTNLP();
      StrengtheningTNLP(const StrengtheningTNLP&);
      StrengtheningTNLP& operator=(const StrengtheningTNLP&);
      //@}

      /** TMINLP (with current bounds) for which the cut it to be
       *  generated */
      const SmartPtr<TMINLP> tminlp_;

      /** Gradient of the (linear) objective function */
      Number* obj_grad_;

      /** Dimension of original problem */
      const Index n_orig_;

      /** Number of constraints in original problem */
      Index m_orig_;

      /** Starting point */
      Number* starting_point_;

      /** Full dimentional x which is used to call the TMINLP
       *  evaluation routines */
      Number* x_full_;

      /** Lower bounds for constraint variables */
      Number* x_l_;

      /** Upper bounds for constraint variables */
      Number* x_u_;

      /** Index of the constraint */
      const Index constr_index_;

      /** Number of variables appearing in the constraint */
      const Index nvar_constr_;

      /** List of variables appearing on the constraints */
      Index* var_indices_;

      /** Flag indicating if the cut has a lower or upper bound */
      bool lower_bound_;

      /** Flag indicating if we TNLP has been solved successfully */
      bool have_final_bound_;

      /** Final strengthened bound */
      Number strengthened_bound_;

      /** space for original gradient if objective function is handled */
      Number* grad_f_;

      /** Auxilliary method for updating the full x variable */
      void update_x_full(const Number *x);
    };

  public:
    /** @name Constructor/Destructor */
    //@{
    /** Constructor.  It is given a TNLP solver to solve the internal
     *  NLPs. */
    CutStrengthener(SmartPtr<TNLPSolver> tnlp_solver,
		    SmartPtr<OptionsList> options);

    /** Destructor */
    virtual ~CutStrengthener();
    //@}

    /** Method for generating and strenghtening all desired cuts */
    bool ComputeCuts(OsiCuts &cs,
		     TMINLP* tminlp,
		     TMINLP2TNLP* problem,
		     const int gindex, CoinPackedVector& cut,
		     double& cut_lb, double& cut_ub,
		     const double g_val, const double g_lb,
		     const double g_ub,
		     int n, const double* x,
		     double infty);

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    CutStrengthener();

    /** Copy Constructor */
    CutStrengthener(const CutStrengthener&);

    /** Overloaded Equals Operator */
    void operator=(const CutStrengthener&);
    //@}

    /** Method for strengthening one cut. */
    bool StrengthenCut(SmartPtr<TMINLP> tminlp /** current TMINLP */,
		       int constr_index /** Index number of the constraint to be strengthened, -1 means objective function */,
		       const CoinPackedVector& row /** Cut to be strengthened */,
		       int n /** Number of variables */,
		       const double* x /** solution from node */,
		       const double* x_l /** Lower bounds for x in which should be valid. */,
		       const double* x_u /** Upper bounds for x in which should be valid. */,
		       double& lb,
		       double& ub);

    /** Method for generating one type of cut (strengthened or disjunctive) */
    bool HandleOneCut(bool is_tight, TMINLP* tminlp,
		      TMINLP2TNLP* problem,
		      const double* minlp_lb,
		      const double* minlp_ub,
		      const int gindex, CoinPackedVector& cut,
		      double& cut_lb, double& cut_ub,
		      int n, const double* x,
		      double infty);

    /** Object for solving the TNLPs */
    SmartPtr<TNLPSolver> tnlp_solver_;

    /** Type of OA cut strengthener */
    int cut_strengthening_type_;
    /** What kind of disjuntion should be done */
    int disjunctive_cut_type_;
    /** verbosity level for OA-related output */
    int oa_log_level_;
  };

} // namespace Ipopt
#endif

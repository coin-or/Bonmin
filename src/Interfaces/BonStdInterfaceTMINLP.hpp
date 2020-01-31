// (C) Copyright Shaurya Sharma
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// S. Sharma
//
// Date :  07.06.2015

#ifndef BONSTDINTERFACETMINLP_HPP
#define BONSTDINTERFACETMINLP_HPP

#include "BonTMINLP.hpp"
#include "BonStdCInterface.h"
#include "IpJournalist.hpp"
#include "IpException.hpp"
#include "IpSmartPtr.hpp"


namespace Bonmin
{
    class StdInterfaceTMINLP : public TMINLP
    {
    public:
        StdInterfaceTMINLP(Ipopt::Index n_var,
                         const Ipopt::Number* x_L, const Ipopt::Number* x_U,
                         Ipopt::Index n_con,
                         const Ipopt::Number* g_L, const Ipopt::Number* g_U,
                         Ipopt::Index nele_jac,
                         Ipopt::Index nele_hess,
                         Ipopt::Index index_style,
                         const Ipopt::Number* start_x,
                         const Ipopt::Number* start_lam,
                         const Ipopt::Number* start_z_L,
                         const Ipopt::Number* start_z_U,
                         Eval_F_CB eval_f,
                         Eval_G_CB eval_g,
                         Eval_Grad_F_CB eval_grad_f,
                         Eval_Jac_G_CB eval_jac_g,
                         Eval_H_CB eval_h,
                         VariableType* var_types,
                         Ipopt::TNLP::LinearityType* var_linearity_types,
                         Ipopt::TNLP::LinearityType* constraint_linearity_types,
                         // BranchingInfo branch,
                         // SosInfo sos,
                         // Intermediate_CB intermediate_cb,
                         Ipopt::Number* x_sol,
                         Ipopt::Number* z_L_sol,
                         Ipopt::Number* z_U_sol,
                         Ipopt::Number* g_sol,
                         Ipopt::Number* lam_sol,
                         Ipopt::Number* obj_sol,
                         UserDataPtr user_data,
                         Ipopt::Number obj_scaling = 1,
                         const Ipopt::Number* x_scaling = NULL,
                         const Ipopt::Number* g_scaling = NULL);



        virtual ~StdInterfaceTMINLP();

        virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
            Ipopt::Index& nnz_h_lag, Ipopt::TNLP::IndexStyleEnum& index_style);

        virtual bool get_scaling_parameters(Ipopt::Number& obj_scaling,
            bool& use_x_scaling, Ipopt::Index n,
            Ipopt::Number* x_scaling,
            bool& use_g_scaling, Ipopt::Index m,
            Ipopt::Number* g_scaling);

        virtual bool get_variables_types(Ipopt::Index n, VariableType* var_types);

        virtual bool get_variables_linearity(Ipopt::Index n, Ipopt::TNLP::LinearityType* var_types);

        virtual bool get_constraints_linearity(Ipopt::Index m, Ipopt::TNLP::LinearityType* const_types);

        virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u);

        virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
            Ipopt::Index m, bool init_lambda,
            Ipopt::Number* lambda);

        virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value);

        virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f);

        virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g);

        virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
          Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
          Ipopt::Index *jCol, Ipopt::Number* values);

        virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
            Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
            bool new_lambda, Ipopt::Index nele_hess,
            Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values);

        virtual void finalize_solution(TMINLP::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x, Ipopt::Number obj_value){}

        void apply_new_x(bool new_x, Ipopt::Index n, const Ipopt::Number* x);

        // TODO
        virtual const SosInfo * sosConstraints() const{return sos_info_;}
        virtual const BranchingInfo* branchingInfo() const{return branch_info_;}

        SosInfo* getSosInfo(){return sos_info_;}
        BranchingInfo* getBranchingInfo(){return branch_info_;}

    private:
        /** Journlist *info_/
        Ipopt::SmartPtr<const Ipopt::Journalist> jnlst_; */

        /** @name Information about the problem */
        //@{
        /** Ipopt::Number of variables */
        const Ipopt::Index n_var_;
        /** Ipopt::Number of constraints */
        const Ipopt::Index n_con_;
        /** Pointer to Ipopt::Number array containing lower bounds for variables */
        const Ipopt::Number* x_L_;
        /** Pointer to Ipopt::Number array containing upper bounds for variables */
        const Ipopt::Number* x_U_;
        /** Pointer to Ipopt::Number array containing lower bounds for constraints */
        const Ipopt::Number* g_L_;
        /** Pointer to Ipopt::Number array containing upper bounds for constraints */
        const Ipopt::Number* g_U_;
        /** Ipopt::Number of non-zero elements in the constraint Jacobian */
        const Ipopt::Index nele_jac_;
        /** Ipopt::Number of non-zero elements in the Hessian */
        const Ipopt::Index nele_hess_;
        /** Starting value of the iRow and jCol parameters for matrices */
        const Ipopt::Index index_style_;
        /** Pointer to Ipopt::Number array containing starting point for variables */
        const Ipopt::Number* start_x_;
        /** Poitner to Ipopt::Number array containing starting values for
        *  constraint multipliers */
        const Ipopt::Number* start_lam_;
        /** Pointer to Ipopt::Number array containing starting values for lower
        *  bound multipliers */
        const Ipopt::Number* start_z_L_;
        /** Pointer to Ipopt::Number array containing starting values for upper
        *  bound multipliers */
        const Ipopt::Number* start_z_U_;
        /** Pointer to callback function evaluating value of objective function */
        Eval_F_CB eval_f_;
        /**  Pointer to callback function evaluating value of constraints */
        Eval_G_CB eval_g_;
        /** Pointer to callback function evaluating gradient of objective
        *  function */
        Eval_Grad_F_CB eval_grad_f_;
        /** Pointer to callback function evaluating Jacobian of constraints */
        Eval_Jac_G_CB eval_jac_g_;

        /** Pointer to callback function evaluating Hessian of Lagrangian */
        Eval_H_CB eval_h_;

        VariableType* var_types_;

        Ipopt::TNLP::LinearityType* var_linearity_types_;

        Ipopt::TNLP::LinearityType* constraint_linearity_types_;

        /** Pointer to intermediate callback function giving control to user */
        Intermediate_CB intermediate_cb_;

        /** Storage of branching priorities information.*/
        BranchingInfo* branch_info_;

        /** Storage of sos constraints */
        SosInfo* sos_info_;

        /** Pointer to user data */
        UserDataPtr user_data_;
        /** Objective scaling factor */
        Ipopt::Number obj_scaling_;
        /** Scaling factors for variables (if not NULL) */
        const Ipopt::Number* x_scaling_;
        /** Scaling factors for constraints (if not NULL) */
        const Ipopt::Number* g_scaling_;
        //@}


        /** A non-const copy of x - this is kept up-to-date in apply_new_x */
        Ipopt::Number* non_const_x_;

        /** Pointers to the user provided vectors for solution */
        Ipopt::Number* x_sol_;
        Ipopt::Number* z_L_sol_;
        Ipopt::Number* z_U_sol_;
        Ipopt::Number* g_sol_;
        Ipopt::Number* lambda_sol_;
        Ipopt::Number* obj_sol_;

        // /** Default Constructor */
        // StdInterfaceTMINLP();

        // /** Copy Constructor */
        // StdInterfaceTMINLP(const StdInterfaceTMINLP&);

        // /** Overloaded Equals Operator */
        // void operator=(const StdInterfaceTMINLP&);
        // //@}
    };
}


#endif // BONSTDINTERFACETMINLP_HPP

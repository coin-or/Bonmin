// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2004
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Carl D. Laird, Carnegie Mellon University,
// Andreas Waechter, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 12/01/2004

#ifndef __IPAMPLTMINLP_HPP__
#define __IPAMPLTMINLP_HPP__

#include "TMINLP.hpp"
#include "IpSmartPtr.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiCuts.hpp"

namespace Ipopt
{


  // Declarations, so that we don't have to include the Ipopt AMPL headers
  class AmplTNLP;
  class AmplSuffixHandler;
  class AmplOptionsList;

  /** Ampl MINLP Interface.
   *  Ampl MINLP Interface, implemented as a TMINLP.
   *  This interface creates a AmplTNLP and also retrieves
   *  the information about the binary and integer variables
   */
  class AmplTMINLP : public TMINLP
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor */
    AmplTMINLP(const SmartPtr<const Journalist>& jnlst,
        const SmartPtr<OptionsList> options,
        char**& argv,
        AmplSuffixHandler* suffix_handler = NULL,
        const std::string& appName = "bonmin",
        std::string* nl_file_content = NULL);

    virtual void Initialize(const SmartPtr<const Journalist>& jnlst,
        const SmartPtr<OptionsList> options,
        char**& argv,
        AmplSuffixHandler* suffix_handler =NULL,
        const std::string& appName = "bonmin",
        std::string* nl_file_content = NULL);

    /** read the branching priorities from ampl suffixes.*/
    void read_priorities();

    /** read the sos constraints from ampl suffixes */
    void read_sos();

    AmplTMINLP();

    virtual AmplTMINLP * createEmpty()
    {
      AmplTMINLP * tminlp = new AmplTMINLP;
      return tminlp;
    }

    /** Default destructor */
    virtual ~AmplTMINLP();
    //@}

    /**@name methods to gather information about the NLP. These
    * methods are overloaded from TMINLP. See TMINLP for their more
    * detailed documentation. */
    //@{
    /** returns dimensions of the nlp. Overloaded from TMINLP */
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
        Index& nnz_h_lag,
        TNLP::IndexStyleEnum& index_style);

    /** returns the vector of variable types */
    virtual bool get_var_types(Index n, VariableType* var_types);

    /** return the variables linearity (linear or not)*/
    virtual bool get_variables_linearity(Index n, Linearity* var_types);
    
    /** return the vector of constraints types*/
    virtual bool get_constraints_types(Index n, Linearity* const_types);
    /** returns bounds of the nlp. Overloaded from TMINLP */
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
        Index m, Number* g_l, Number* g_u);

    /** provides a starting point for the nlp variables. Overloaded
    from TMINLP */
    virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
        Index m, bool init_lambda, Number* lambda);

    /** evaluates the objective value for the nlp. Overloaded from TMINLP */
    virtual bool eval_f(Index n, const Number* x, bool new_x,
        Number& obj_value);

    /** evaluates the gradient of the objective for the
    nlp. Overloaded from TMINLP */
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x,
        Number* grad_f);

    /** evaluates the constraint residuals for the nlp. Overloaded from TMINLP */
    virtual bool eval_g(Index n, const Number* x, bool new_x,
        Index m, Number* g);

    /** specifies the jacobian structure (if values is NULL) and
     *  evaluates the jacobian values (if values is not NULL) for the
     *  nlp. Overloaded from TMINLP */
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
        Index m, Index nele_jac, Index* iRow,
        Index *jCol, Number* values);

    /** specifies the structure of the hessian of the lagrangian (if
     *  values is NULL) and evaluates the values (if values is not
     *  NULL). Overloaded from TMINLP */
    virtual bool eval_h(Index n, const Number* x, bool new_x,
        Number obj_factor, Index m, const Number* lambda,
        bool new_lambda, Index nele_hess, Index* iRow,
        Index* jCol, Number* values);
    //@}

    /** @name Solution Methods */
    //@{
    virtual void finalize_solution(SolverReturn status,
        Index n, const Number* x, Number obj_value) const ;

    void write_solution(const std::string & message, const Number *x_sol) const;
    //@}

    /** Write the solution file.  This is a wrapper for AMPL's
     *  write_sol.  TODO Maybe this should be at a different place, or
     *  collect the numbers itself? */
    void write_solution_file(const std::string& message, const Number * x) const;
    //@}

   
    virtual const BranchingInfo * branchingInfo() const
    {
      return &branch_;
    } 

    virtual const SosInfo * sosConstraints() const
    {
      return &sos_;
    }

    /** @name User callbacks */
    //@{
    /** Additional application specific options.*/
    virtual void fillApplicationOptions(AmplOptionsList* amplOptList)
    {}
    //@}


    /** This methods gives the linear part of the objective function */
    virtual void getLinearPartOfObjective(double * obj);
    /** Journlist */
    SmartPtr<const Journalist> jnlst_;

    /** Sign of the objective fn (1 for min, -1 for max) */
    double obj_sign_;

    /**@name Problem Size Data*/
    //@{
    Index nz_h_full_; // number of nonzeros in the full_x hessian
    /* the rest of the problem size data is available easily through the ampl variables */
    //@}

    /**@name Internal copies of data */
    //@{
    /** A non-const copy of x - this is kept up-to-date in apply_new_x */
    Number* non_const_x_;

    /** Solution Vectors */
    Number* x_sol_;
    Number* z_L_sol_;
    Number* z_U_sol_;
    Number* g_sol_;
    Number* lambda_sol_;
    Number obj_sol_;
    //@}

    /**@name Flags to track internal state */
    //@{
    /** true when the objective value has been calculated with the
     *  current x, set to false in apply_new_x, and set to true in
     *  internal_objval */
    bool objval_called_with_current_x_;
    /** true when the constraint values have been calculated with the
     *  current x, set to false in apply_new_x, and set to true in
     *  internal_conval */
    bool conval_called_with_current_x_;
    //@}


    /** Make the objective call to ampl */
    bool internal_objval(Number& obj_val);

    /** Make the constraint call to ampl*/
    bool internal_conval(Index m, Number* g=NULL);

    /** Internal function to update the internal and ampl state if the
    x value changes */
    void apply_new_x(bool new_x, Index n, const Number* x);

    /** Method to add the extra option for ampl */
    void fillAmplOptionList(AmplOptionsList* amplOptList);

    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{

    /** Copy Constructor */
    AmplTMINLP(const AmplTMINLP&);

    /** Overloaded Equals Operator */
    void operator=(const AmplTMINLP&);
    //@}

    /** pointer to the internal AmplTNLP */
    AmplTNLP* ampl_tnlp_;

    /** Storage of branching priorities information.*/
    BranchingInfo branch_;
    /** Storage of sos constraints */
    SosInfo sos_;
    /** Store a suffix handler (since Ipopt does not want to give access to his >:( ).*/
    SmartPtr<AmplSuffixHandler> suffix_handler_;
  };
} // namespace Ipopt

#endif

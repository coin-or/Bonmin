// (C) Copyright International Business Machines Corporation and
// Carnegie Mellon University 2004, 2007
//
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Carl D. Laird, Carnegie Mellon University,
// Andreas Waechter, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 12/01/2004

#ifndef __IPAMPLTMINLP_HPP__
#define __IPAMPLTMINLP_HPP__

#include "BonTMINLP.hpp"
#include "IpSmartPtr.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiCuts.hpp"
#include "BonRegisteredOptions.hpp"
#include "BonTypes.hpp"

/* non Ipopt forward declaration */
struct ASL_pfgh;
struct SufDecl;
struct SufDesc;


// Declarations, so that we don't have to include the Ipopt AMPL headers
namespace Ipopt
{
  class AmplSuffixHandler;
  class AmplOptionsList;
  class AmplTNLP;
}

namespace Bonmin
{

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
    AmplTMINLP(const Ipopt::SmartPtr<const Ipopt::Journalist>& jnlst,
        const Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
        const Ipopt::SmartPtr<Ipopt::OptionsList> options,
        char**& argv,
        Ipopt::AmplSuffixHandler* suffix_handler = NULL,
        const std::string& appName = "bonmin",
        std::string* nl_file_content = NULL);

    virtual void Initialize(const Ipopt::SmartPtr<const Ipopt::Journalist>& jnlst,
        const Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
        const Ipopt::SmartPtr<Ipopt::OptionsList> options,
        char**& argv,
        Ipopt::AmplSuffixHandler* suffix_handler =NULL,
        const std::string& appName = "bonmin",
        std::string* nl_file_content = NULL);

    /** read the branching priorities from ampl suffixes.*/
    void read_priorities();

    /** read the sos constraints from ampl suffixes */
    void read_sos();

    /** Read suffixes which indicate which constraints are convex.*/
    void read_convexities();

    /** Read suffixes used to apply perspective in OA to some of the constraints.*/
    void read_onoff();

    /** Read suffixes on objective functions for upper bounding*/
    void read_obj_suffixes();

    /** Default constructor.*/
    AmplTMINLP();

    virtual AmplTMINLP * createEmpty()
    {
      AmplTMINLP * tminlp = new AmplTMINLP;
      return tminlp;
    }

    /** destructor */
    virtual ~AmplTMINLP();
    //@}

    /** Return the ampl solver object (ASL*) */
    const ASL_pfgh* AmplSolverObject() const;


    /**@name methods to gather information about the NLP. These
    * methods are overloaded from TMINLP. See TMINLP for their more
    * detailed documentation. */
    //@{
    /** returns dimensions of the nlp. Overloaded from TMINLP */
    virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag,
        Ipopt::TNLP::IndexStyleEnum& index_style);

    /** returns the vector of variable types */
    virtual bool get_variables_types(Ipopt::Index n, VariableType* var_types);

    /** return the variables linearity (linear or not)*/
    virtual bool get_variables_linearity(Ipopt::Index n, Ipopt::TNLP::LinearityType * var_types);

    /** Returns the constraint linearity.
     * array should be alocated with length at least n.*/
    virtual bool get_constraints_linearity(Ipopt::Index m,
        Ipopt::TNLP::LinearityType* const_types);

    /** returns bounds of the nlp. Overloaded from TMINLP */
    virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
        Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u);

    /** provides a starting point for the nlp variables. Overloaded
    from TMINLP */
    virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
        bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
        Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda);

    /** evaluates the objective value for the nlp. Overloaded from TMINLP */
    virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
        Ipopt::Number& obj_value);

    /** evaluates the gradient of the objective for the
    nlp. Overloaded from TMINLP */
    virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
        Ipopt::Number* grad_f);

    /** evaluates the constraint residuals for the nlp. Overloaded from TMINLP */
    virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
        Ipopt::Index m, Ipopt::Number* g);

    /** specifies the jacobian structure (if values is NULL) and
     *  evaluates the jacobian values (if values is not NULL) for the
     *  nlp. Overloaded from TMINLP */
    virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
        Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
        Ipopt::Index *jCol, Ipopt::Number* values);

    /** specifies the structure of the hessian of the lagrangian (if
     *  values is NULL) and evaluates the values (if values is not
     *  NULL). Overloaded from TMINLP */
    virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
        Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
        bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
        Ipopt::Index* jCol, Ipopt::Number* values);

    /** compute the value of a single constraint */
    virtual bool eval_gi(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
        Ipopt::Index i, Ipopt::Number& gi);
    /** compute the structure or values of the gradient for one
    constraint */
    virtual bool eval_grad_gi(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
        Ipopt::Index i, Ipopt::Index& nele_grad_gi, Ipopt::Index* jCol,
        Ipopt::Number* values);
    //@}

    /** @name Solution Methods */
    //@{
    /** Called after optimizing to return results to ampl.
      * Status code is put into solve_result_num according to the table below.
      * <table>
      * <tr> <td> <b> <center> Code </center> </b> </td> <td> <b> <center> Status </center> </b> </td> </tr>
      * <tr> <td> 3 </td> <td> Integer optimal </td> </tr>
      * <tr> <td> 220 </td> <td> problem is proven infeasible. </td> </tr>
      * <tr> <td> 421 </td> <td> limit reached or user interrupt with integer feasible solution found. </td> </tr>
      * <tr> <td> 410 </td> <td> limit reached or user interrupt without any integer feasible solution. </td> </tr>
      * <tr> <td> 500 </td> <td> error. </td> </tr>
      * <caption> Status codes for optimization. </caption>
      * </table>
      *    */
    virtual void finalize_solution(TMINLP::SolverReturn status,
        Ipopt::Index n, const Ipopt::Number* x, Ipopt::Number obj_value);

    /** Write the solution using ampl's write_sol (called by finalize_solution).*/
    void write_solution(const std::string & message, const Ipopt::Number *x_sol);
    //@}

    //@}


    virtual const BranchingInfo * branchingInfo() const
    {
      return &branch_;
    }

    virtual const SosInfo * sosConstraints() const
    {
      return &sos_;
    }

    virtual const PerturbInfo* perturbInfo() const
    {
      return &perturb_info_;
    }

    /** @name User callbacks */
    //@{
    /** Additional application specific options.*/
    virtual void fillApplicationOptions(Ipopt::AmplOptionsList* amplOptList)
    {}
    //@}


    /** This methods gives the linear part of the objective function */
    virtual void getLinearPartOfObjective(double * obj);


    /** Do we have an alternate objective for upper bounding?*/
    virtual bool hasUpperBoundingObjective()
    {
      return upperBoundingObj_ != -1;
    }

    /** This method to returns the value of an alternative objective function for
      upper bounding (if one has been declared by using the prefix UBObj).*/
    virtual bool eval_upper_bound_f(Ipopt::Index n, const Ipopt::Number* x,
        Ipopt::Number& obj_value);

    /** Get accest to constraint convexities.*/
    virtual bool get_constraint_convexities(int m, TMINLP::Convexity * constraints_convexities)const
    {
      if (constraintsConvexities_ != NULL) {
        CoinCopyN(constraintsConvexities_, m, constraints_convexities);
      }
      else {
        CoinFillN(constraints_convexities, m, TMINLP::Convex);
      }
      return true;
    }
    /** Get dimension information on nonconvex constraints.*/
    virtual bool get_number_nonconvex(int & number_non_conv, int & number_concave) const
    {
      number_non_conv = numberNonConvex_;
      number_concave = numberSimpleConcave_;
      return true;
    }
    /** Get array describing the constraints marked nonconvex in the model.*/
    virtual bool get_constraint_convexities(int number_non_conv, MarkedNonConvex * non_convexes) const
    {
      assert(number_non_conv == numberNonConvex_);
      CoinCopyN( nonConvexConstraintsAndRelaxations_, number_non_conv, non_convexes);
      return true;
    }
    /** Fill array containing indices of simple concave constraints.*/
    virtual bool get_simple_concave_constraints(int number_concave, SimpleConcaveConstraint * simple_concave) const
    {
      assert(number_concave == numberSimpleConcave_);
      CoinCopyN(simpleConcaves_, numberSimpleConcave_, simple_concave);
      return true;
    }

    /** Say if problem has a linear objective (for OA) */
    virtual bool hasLinearObjective()
    {
      return hasLinearObjective_;
    }

  /** Access array describing onoff constraint.*/
  virtual const int * get_const_xtra_id() const{
    return c_extra_id_();
  }
  private:
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
     /** Name of application.*/
    std::string appName_;

    /** Index of the objective to use for upper bounding*/
    int upperBoundingObj_;
    /** pointer to the internal AmplTNLP */
    Ipopt::AmplTNLP* ampl_tnlp_;
    /** Journalist */
    Ipopt::SmartPtr<const Ipopt::Journalist> jnlst_;

    /** Storage of branching priorities information.*/
    BranchingInfo branch_;
    /** Storage of sos constraints */
    SosInfo sos_;
    /** Storage for perturbation radii */
    PerturbInfo perturb_info_;
    /** Store a suffix handler.*/
    Ipopt::SmartPtr<Ipopt::AmplSuffixHandler> suffix_handler_;

    /** Store constraints types.*/
    TMINLP::Convexity * constraintsConvexities_;

    /** Store onoff information.*/
    vector<int> c_extra_id_; 

    /** Ipopt::Number of nonConvex constraints.*/
    int numberNonConvex_;
    /** Store marked non-convex constraints and their relaxations.*/
    MarkedNonConvex * nonConvexConstraintsAndRelaxations_;
    /** Ipopt::Number of simpleConcave constraints.*/
    int numberSimpleConcave_;
    /** Store simple concave constraints descriptions.*/
    SimpleConcaveConstraint * simpleConcaves_;

    /** Flag to indicate if objective function is linear */
    bool hasLinearObjective_;

    /** Flag to say if AMPL solution file should be written.*/
    int writeAmplSolFile_;
  };
} // namespace Ipopt

#endif


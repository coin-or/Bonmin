// (C) Copyright International Business Machines Corporation and
// Carnegie Mellon University 2004, 2007
//
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

#include "BonTMINLP.hpp"
#include "IpSmartPtr.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiCuts.hpp"
#include "BonRegisteredOptions.hpp"

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
    AmplTMINLP(const SmartPtr<const Journalist>& jnlst,
        const SmartPtr<Bonmin::RegisteredOptions> roptions,
        const SmartPtr<OptionsList> options,
        char**& argv,
        AmplSuffixHandler* suffix_handler = NULL,
        const std::string& appName = "bonmin",
        std::string* nl_file_content = NULL);

    virtual void Initialize(const SmartPtr<const Journalist>& jnlst,
        const SmartPtr<Bonmin::RegisteredOptions> roptions,
        const SmartPtr<OptionsList> options,
        char**& argv,
        AmplSuffixHandler* suffix_handler =NULL,
        const std::string& appName = "bonmin",
        std::string* nl_file_content = NULL);

    /** read the branching priorities from ampl suffixes.*/
    void read_priorities();

    /** read the sos constraints from ampl suffixes */
    void read_sos();

    /** Read suffixes which indicate which constraints are convex.*/
    void read_convexities();

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
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
        Index& nnz_h_lag,
        TNLP::IndexStyleEnum& index_style);

    /** returns the vector of variable types */
    virtual bool get_variables_types(Index n, VariableType* var_types);

    /** return the variables linearity (linear or not)*/
    virtual bool get_variables_linearity(Index n, Ipopt::TNLP::LinearityType * var_types);

    /** Returns the constraint linearity.
     * array should be alocated with length at least n.*/
    virtual bool get_constraints_linearity(Index m,
        Ipopt::TNLP::LinearityType* const_types);

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

    /** compute the value of a single constraint */
    virtual bool eval_gi(Index n, const Number* x, bool new_x,
        Index i, Number& gi);
    /** compute the structure or values of the gradient for one
    constraint */
    virtual bool eval_grad_gi(Index n, const Number* x, bool new_x,
        Index i, Index& nele_grad_gi, Index* jCol,
        Number* values);
    //@}

    /** @name Solution Methods */
    //@{
    virtual void finalize_solution(TMINLP::SolverReturn status,
        Index n, const Number* x, Number obj_value);

    void write_solution(const std::string & message, const Number *x_sol);
    //@}

    /** Write the solution file.  This is a wrapper for AMPL's
     *  write_sol.  TODO Maybe this should be at a different place, or
     *  collect the numbers itself? */
    void write_solution_file(const std::string& message) const;
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
    virtual void fillApplicationOptions(AmplOptionsList* amplOptList)
    {}
    //@}


    /** This methods gives the linear part of the objective function */
    virtual void getLinearPartOfObjective(double * obj);


    /** Method to add the extra option for ampl */
    void fillAmplOptionList(AmplOptionsList* amplOptList);

    /** Do we have an alternate objective for upper bounding?*/
    virtual bool hasUpperBoundingObjective()
    {
      return upperBoundingObj_ != -1;
    }

    /** This method to returns the value of an alternative objective function for
      upper bounding (if one has been declared by using the prefix UBObj).*/
    virtual bool eval_upper_bound_f(Index n, const Number* x,
        Number& obj_value);

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

    /** Index of the objective to use for upper bounding*/
    int upperBoundingObj_;
    /** pointer to the internal AmplTNLP */
    AmplTNLP* ampl_tnlp_;
    /** Journalist */
    SmartPtr<const Journalist> jnlst_;

    /** Storage of branching priorities information.*/
    BranchingInfo branch_;
    /** Storage of sos constraints */
    SosInfo sos_;
    /** Storage for perturbation radii */
    PerturbInfo perturb_info_;
    /** Store a suffix handler.*/
    SmartPtr<AmplSuffixHandler> suffix_handler_;

    /** Store constraints types.*/
    TMINLP::Convexity * constraintsConvexities_;

    /** Number of nonConvex constraints.*/
    int numberNonConvex_;
    /** Store marked non-convex constraints and their relaxations.*/
    MarkedNonConvex * nonConvexConstraintsAndRelaxations_;
    /** Number of simpleConcave constraints.*/
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


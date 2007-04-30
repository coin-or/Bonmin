// (C) Copyright International Business Machines Corporation and
// Carnegie Mellon University 2004, 2007
//
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
// Carl D. Laird, Carnegie Mellon University,
// Andreas Waechter, International Business Machines Corporation
//
// Date : 12/01/2004

#ifndef __TMINLP_HPP__
#define __TMINLP_HPP__

#include "IpUtils.hpp"
#include "IpReferenced.hpp"
#include "IpException.hpp"
#include "IpAlgTypes.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiCuts.hpp"
#include "IpTNLP.hpp"
#include "CoinError.hpp"

#include "CoinHelperFunctions.hpp"
using namespace Ipopt;
namespace Bonmin
{
  DECLARE_STD_EXCEPTION(TMINLP_INVALID);
  DECLARE_STD_EXCEPTION(TMINLP_INVALID_VARIABLE_BOUNDS);

  /** Base class for all MINLPs that use a standard triplet matrix form
   *  and dense vectors.
   *  The class TMINLP2TNLP allows the caller to produce a viable TNLP
   *  from the MINLP (by relaxing binary and/or integers, or by
   *  fixing them), which can then be solved by Ipopt.
   *
   *  This interface presents the problem form:
   *  \f[
   *  \begin{array}{rl}
   *     &min f(x)\\
   *
   *     \mbox{s.t.}&\\
   *      &   g^L <= g(x) <= g^U\\
   *
   *       &   x^L <=  x   <= x^U\\
   *   \end{array}
   *  \f]
   *  Where each x_i is either a continuous, binary, or integer variable.
   *  If x_i is binary, the bounds [xL,xU] are assumed to be [0,1].
   *  In order to specify an equality constraint, set gL_i = gU_i =
   *  rhs.  The value that indicates "infinity" for the bounds
   *  (i.e. the variable or constraint has no lower bound (-infinity)
   *  or upper bound (+infinity)) is set through the option
   *  nlp_lower_bound_inf and nlp_upper_bound_inf.  To indicate that a
   *  variable has no upper or lower bound, set the bound to
   *  -ipopt_inf or +ipopt_inf respectively
   */
  class TMINLP : public Ipopt::ReferencedObject
  {
  public:
    friend class TMINLP2TNLP;
    /** Return statuses of algorithm.*/
    enum SolverReturn{
      SUCCESS,
      INFEASIBLE,
      LIMIT_EXCEEDED,
      MINLP_ERROR};
    /** Class to store sos constraints for model */
    struct SosInfo
    {
      /** Number of SOS constraints.*/
      int num;
      /** Type of sos. At present Only type '1' SOS are supported by Cbc*/
      char * types;
      /** priorities of sos constraints.*/
      int * priorities;
      
      /** \name Sparse storage of the elements of the SOS constraints.*/
      /** @{ */
      /** Total number of non zeroes in SOS constraints.*/
      int numNz;
      /** For 0 <= i < nums, start[i] gives the indice of indices and weights arrays at which the description of constraints i begins..*/ 
      int * starts;
      /** indices of elements belonging to the SOS.*/
      int * indices;
      /** weights of the elements of the SOS.*/
      double * weights;
      /** @} */
      /** default constructor. */
      SosInfo();
      /** Copy constructor.*/
      SosInfo(const SosInfo & source);
      

      /** destructor*/
      ~SosInfo()
      {
        gutsOfDestructor();
      }


      /** Reset information */
      void gutsOfDestructor();

    };

    /** Stores branching priorities information. */
    struct BranchingInfo
    {
      /**number of variables*/
      int size;
      /** User set priorities on variables. */
      int * priorities;
      /** User set preferered branching direction. */
      int * branchingDirections;
      /** User set up pseudo costs.*/
      double * upPsCosts;
      /** User set down pseudo costs.*/
      double * downPsCosts;
      BranchingInfo():
      size(0),
      priorities(NULL),
      branchingDirections(NULL),
      upPsCosts(NULL),
      downPsCosts(NULL)
      {}
      BranchingInfo(const BranchingInfo &other)
      {
        gutsOfDestructor();
        size = other.size;
        priorities = CoinCopyOfArray(other.priorities, size);
        branchingDirections = CoinCopyOfArray(other.branchingDirections, size);
        upPsCosts = CoinCopyOfArray(other.upPsCosts, size);
        downPsCosts = CoinCopyOfArray(other.downPsCosts, size);
      }
      void gutsOfDestructor()
      {
      if (priorities != NULL) delete [] priorities;
      priorities = NULL;
      if (branchingDirections != NULL) delete [] branchingDirections;  
      branchingDirections = NULL;
      if (upPsCosts != NULL) delete [] upPsCosts;
      upPsCosts = NULL;
      if (downPsCosts != NULL) delete [] downPsCosts;
      downPsCosts = NULL;
      }
      ~BranchingInfo()
      {
	gutsOfDestructor();
      }
    };

    /** Class to store perturbation radii for variables in the model */
    class PerturbInfo
    {
    public:
      /** default constructor. */
      PerturbInfo() :
	perturb_radius_(NULL)
      {}

      /** destructor*/
      ~PerturbInfo()
      {
        delete [] perturb_radius_;
      }

      /** Method for setting the perturbation radii. */
      void SetPerturbationArray(Index numvars, const double* perturb_radius);

      /** Method for getting the array for the perturbation radii in
       *  order to use the values. */
      const double* GetPerturbationArray() const {
	return perturb_radius_;
      }

    private:
      /** Copy constructor.*/
      PerturbInfo(const PerturbInfo & source);

      /** Perturbation radii for all variables.  A negative value
       *  means that the radius has not been given. If the pointer is
       *  NULL, then no variables have been assigned a perturbation
       *  radius. */
      double* perturb_radius_;
    };

    /** Type of the variables.*/
    enum VariableType
    {
      CONTINUOUS,
      BINARY,
      INTEGER
    };

    /**@name Constructors/Destructors */
    //@{
    TMINLP();

    /** Default destructor */
    virtual ~TMINLP()
    {}
    //@}

    /**@name methods to gather information about the MINLP */
    //@{
    /** overload this method to return the number of variables
     *  and constraints, and the number of non-zeros in the jacobian and
     *  the hessian. */
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
        Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)=0;

    /** overload this method to return scaling parameters. This is
     *  only called if the options are set to retrieve user scaling.
     *  There, use_x_scaling (or use_g_scaling) should get set to true
     *  only if the variables (or constraints) are to be scaled.  This
     *  method should return true only if the scaling parameters could
     *  be provided.
     */
    virtual bool get_scaling_parameters(Number& obj_scaling,
                                        bool& use_x_scaling, Index n,
                                        Number* x_scaling,
                                        bool& use_g_scaling, Index m,
                                        Number* g_scaling)
    {
      return false;
    }


    /** overload this method to set the variable type. The var_types
     *  array will be allocated with length n. */
    virtual bool get_variables_types(Index n, VariableType* var_types)=0;

    /** overload this method to return the constraint linearity.
     * array should be alocated with length at least n.*/
    virtual bool get_constraints_linearity(Index m, 
					   Ipopt::TNLP::LinearityType* const_types) = 0;

    /** overload this method to return the information about the bound
     *  on the variables and constraints. The value that indicates
     *  that a bound does not exist is specified in the parameters
     *  nlp_lower_bound_inf and nlp_upper_bound_inf.  By default,
     *  nlp_lower_bound_inf is -1e19 and nlp_upper_bound_inf is
     *  1e19.
     *  An exception will be thrown if x_l and x_u are not 0,1 for binary variables
     */
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
        Index m, Number* g_l, Number* g_u)=0;

    /** overload this method to return the starting point. The bools
     *  init_x and init_lambda are both inputs and outputs. As inputs,
     *  they indicate whether or not the algorithm wants you to
     *  initialize x and lambda respectively. If, for some reason, the
     *  algorithm wants you to initialize these and you cannot, set
     *  the respective bool to false.
     */
    virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
        Index m, bool init_lambda,
        Number* lambda)=0;

    /** overload this method to return the value of the objective function */
    virtual bool eval_f(Index n, const Number* x, bool new_x,
        Number& obj_value)=0;

    /** overload this method to return the vector of the gradient of
     *  the objective w.r.t. x */
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x,
        Number* grad_f)=0;

    /** overload this method to return the vector of constraint values */
    virtual bool eval_g(Index n, const Number* x, bool new_x,
        Index m, Number* g)=0;

    /** overload this method to return the jacobian of the
     *  constraints. The vectors iRow and jCol only need to be set
     *  once. The first call is used to set the structure only (iRow
     *  and jCol will be non-NULL, and values will be NULL) For
     *  subsequent calls, iRow and jCol will be NULL. */
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
        Index m, Index nele_jac, Index* iRow,
        Index *jCol, Number* values)=0;

    /** overload this method to return the hessian of the
     *  lagrangian. The vectors iRow and jCol only need to be set once
     *  (during the first call). The first call is used to set the
     *  structure only (iRow and jCol will be non-NULL, and values
     *  will be NULL) For subsequent calls, iRow and jCol will be
     *  NULL. This matrix is symmetric - specify the lower diagonal
     *  only */
    virtual bool eval_h(Index n, const Number* x, bool new_x,
        Number obj_factor, Index m, const Number* lambda,
        bool new_lambda, Index nele_hess,
        Index* iRow, Index* jCol, Number* values)=0;
    /** Compute the value of a single constraint. The constraint
     *  number is i (starting counting from 0. */
    virtual bool eval_gi(Index n, const Number* x, bool new_x,
			 Index i, Number& gi)
    {
      std::cerr << "Method eval_gi not overloaded from TMINLP\n";
      return false;
    }
    /** Compute the structure or values of the gradient for one
     *  constraint. The constraint * number is i (starting counting
     *  from 0.  Other things are like with eval_jac_g. */
    virtual bool eval_grad_gi(Index n, const Number* x, bool new_x,
			      Index i, Index& nele_grad_gi, Index* jCol,
			      Number* values)
    {
      std::cerr << "Method eval_grad_gi not overloaded from TMINLP\n";
      return false;
    }
    //@}

    /** @name Solution Methods */
    //@{
    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(TMINLP::SolverReturn status,
                                   Index n, const Number* x, Number obj_value) =0;
    //@}
    
    virtual const BranchingInfo * branchingInfo() const = 0;

    /** Add some linear cuts to the problem formulation */
    void addCuts(int numberCuts, const OsiRowCut ** cuts);
  
    /** Remove some cuts to the formulation */
    void removeCuts(int number ,const int * toRemove);

    /** remove the last number cuts.*/
    void removeLastCuts(int number);

    virtual const SosInfo * sosConstraints() const = 0;

    virtual const PerturbInfo* perturbInfo() const
    {
      return NULL;
    }

    /** Say if has a specific function to compute upper bounds*/
    virtual bool hasUpperBoundingObjective(){
      return false;}
    
    /** overload this method to return the value of an alternative objective function for
      upper bounding (to use it hasUpperBoundingObjective should return true).*/
    virtual bool eval_upper_bound_f(Index n, const Number* x,
                                    Number& obj_value){
    return 0;}
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
    //TMINLP();

    /** Copy Constructor */
    TMINLP(const TMINLP&);

    /** Overloaded Equals Operator */
    void operator=(const TMINLP&);
    //@}

  /** resize arrays for linear cuts */
  void resizeLinearCuts(int newNumberCuts, int newNnz);
  /** columnindices of linear cuts. */
   int * jCol_;
  /** rows indices of linear cuts. */
   int * iRow_;
  /** elements of linear cuts.*/
  double * elems_;
  /** lower bounds for linear cuts. */
  double * lower_;
  /** upper bounds for linear cuts. */
  double * upper_;
  /** number of linear cuts.*/
  int nLinearCuts_;
  /** number of non-zeroes in linear cuts*/
  int linearCutsNnz_;
  /** storage size for linear cuts number of cuts.*/
  int linearCutsCapacity_;
  /** storage size for linear cuts number of nnz.*/
  int linearCutsNnzCapacity_;
  };

} // namespace Ipopt

#endif


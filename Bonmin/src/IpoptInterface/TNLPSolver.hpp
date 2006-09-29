// (C) Copyright International Business Machines (IBM) 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, IBM
//
// Date : 26/09/2006


#ifndef TNLPSolver_H
#define TNLPSolver_H
#include "IpTNLP.hpp"

//Some declarations
#include "IpOptionsList.hpp"
#include "IpRegOptions.hpp"

namespace Bonmin  {
/** This is a generic class for calling an NLP solver to solve a TNLP.
    A TNLPSolver is able to solve and resolve a problem, it has some options (stored
    with Ipopt OptionList structure and registeredOptions) it produces some statistics (in SolveStatisctics and sometimes some errorCodes.
*/
class TNLPSolver: public Ipopt::ReferencedObject{
 public:

  enum ReturnStatus /** Standard return statuses for a solver*/{
    iterationLimit/** Solver reached iteration limit. */,
    computationError /** Some error was made in the computations. */,
    illDefinedProblem /** The solver finds that the problem is not well defined. */,
    illegalOption /** An option is not valid. */,
    externalException /** Some unrecovered exception occured in an external tool used by the solver. */,
    exception /** Some unrocevered exception */,
    solvedOptimal /** Problem solved to an optimal solution.*/,
    solvedOptimalTol /** Problem solved to "acceptable level of tolerance. */,
    provenInfeasible /** Infeasibility Proven. */,
    unbounded /** Problem is unbounded.*/
  };
  /// Constructor
   TNLPSolver();

  ///virtual constructor
  virtual TNLPSolver * createNew() = 0;
  

   /// Virtual destructor
   virtual ~TNLPSolver();

   /** Initialize the TNLPSolver (read options from params_file)
   */
   virtual void Initialize(std::string params_file) = 0;

   /** Initialize the TNLPSolver (read options from istream is)
   */
   virtual void Initialize(std::istream& is) = 0;

   /** @name Solve methods */
   //@{
   /// Solves a problem expresses as a TNLP 
   virtual ReturnStatus OptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> & tnlp) = 0;

   /// Resolves a problem expresses as a TNLP 
   virtual ReturnStatus ReOptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> & tnlp) = 0;
   //@}

   ///Get a pointer to RegisteredOptions (generally used to add new ones)
   virtual Ipopt::SmartPtr<Ipopt::RegisteredOptions> RegOptions() = 0;

   /// Get the options (for getting theur values).
   virtual Ipopt::SmartPtr<const Ipopt::OptionsList> Options() const = 0;

   /// Get the options (for getting and setting their values).
   virtual Ipopt::SmartPtr<Ipopt::OptionsList> Options() = 0;

   /// Get the CpuTime of the last optimization.
   virtual double CPUTime() = 0;

   /// Get the iteration count of the last optimization.
   virtual int IterationCount() = 0;

protected:
  bool zeroDimension(const Ipopt::SmartPtr<Ipopt::TNLP> &tnlp, 
		     ReturnStatus &optimization_status);
   private:
   /// There is no copy constructor for this class
   TNLPSolver(TNLPSolver &other); 
};
}
#endif



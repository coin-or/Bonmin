// (C) Copyright International Business Machines (IBM) 2005
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, IBM
//
// Date : 26/09/2006

#ifndef IpoptSolver_HPP
#define IpoptSolver_HPP
#include "BonTNLPSolver.hpp"
#include "IpIpoptApplication.hpp"


namespace Bonmin{
class IpoptSolver: public TNLPSolver {
public:
class UnsolvedIpoptError: public TNLPSolver::UnsolvedError
{
 public:
  UnsolvedIpoptError(int errorNum = 10000):
  TNLPSolver::UnsolvedError(errorNum)
  {}
  virtual const std::string& errorName() const;
 
  virtual const std::string& solverName() const; 
  virtual ~UnsolvedIpoptError(){}
  private:
  static std::string errorNames [17];
  static std::string solverName_;
};

  virtual UnsolvedError * newUnsolvedError(int num){
    return new UnsolvedIpoptError(num);}



  /// Constructor
  IpoptSolver():
    app_(),
    warmStartStrategy_(1)
  {}

  ///virtual constructor
  virtual TNLPSolver * createNew();


  /// Virtual destructor
  virtual ~IpoptSolver();

  /** Initialize the TNLPSolver (read options from params_file)
  */
  virtual void Initialize(std::string params_file);

  /** Initialize the TNLPSolver (read options from istream is)
  */
  virtual void Initialize(std::istream& is);

  /** @name Solve methods */
  //@{
  /// Solves a problem expresses as a TNLP 
  virtual TNLPSolver::ReturnStatus OptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> & tnlp);

  /// Resolves a problem expresses as a TNLP 
  virtual TNLPSolver::ReturnStatus ReOptimizeTNLP(const Ipopt::SmartPtr<Ipopt::TNLP> & tnlp);

  /// Set the warm start in the solver
  virtual bool setWarmStart(const CoinWarmStart * warm,
                           Ipopt::SmartPtr<TMINLP2TNLP> tnlp);

  /// Get the warm start form the solver
  virtual CoinWarmStart * getWarmStart(Ipopt::SmartPtr<Bonmin::TMINLP2TNLP> tnlp) const;

  virtual CoinWarmStart * getEmptyWarmStart() const;

  /// Enable the warm start options in the solver
  virtual void enableWarmStart();

  /// Disable the warm start options in the solver
  virtual void disableWarmStart();

  //@}

  ///Get a pointer to RegisteredOptions (generally used to add new ones)
  virtual Ipopt::SmartPtr<Ipopt::RegisteredOptions> RegOptions();

  /// Get the options (for getting their values).
  virtual Ipopt::SmartPtr<const Ipopt::OptionsList> Options() const;

   /// Get the options (for getting and setting their values).
   virtual Ipopt::SmartPtr<Ipopt::OptionsList> Options();

   /// Get the CpuTime of the last optimization.
   virtual double CPUTime();

   /// Get the iteration count of the last optimization.
   virtual int IterationCount();

  /// turn off all output from the solver 
  virtual void turnOffOutput();
  /// turn on all output from the solver
  virtual void turnOnOutput();

  /// Return status of last optimization
  Ipopt::ApplicationReturnStatus getOptStatus() const
  {
    return optimizationStatus_;
  }

  Ipopt::IpoptApplication& getIpoptApp(){
    return app_;
  }
  private:
  /** Set default Ipopt parameters for use in a MINLP */
  void setMinlpDefaults(Ipopt::SmartPtr< Ipopt::OptionsList> Options);

  /** get Bonmin return status from Ipopt one. */
  TNLPSolver::ReturnStatus solverReturnStatus(Ipopt::ApplicationReturnStatus optimization_status) const;
  /** Ipopt application */
  Ipopt::IpoptApplication app_;
  /** return status of last optimization.*/
  Ipopt::ApplicationReturnStatus optimizationStatus_;
    //@}

  /** Warm start strategy :
  <ol>
  <li> no warm start,</li>
  <li> simple warm start (optimal point),</li>
  <li> more elaborate strategies (interior point...).</li>
  </ol>
  */
  int warmStartStrategy_;
};
}
#endif


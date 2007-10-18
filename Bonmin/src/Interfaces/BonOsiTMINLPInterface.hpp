// (C) Copyright International Business Machines Corporation, Carnegie Mellon University 2004, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
// Carl D. Laird, Carnegie Mellon University,
// Andreas Waechter, International Business Machines Corporation
//
// Date : 12/01/2004


#ifndef OsiTMINLPInterface_H
#define OsiTMINLPInterface_H

#define INT_BIAS 0e-8

#include <string>
#include <iostream>

#include "OsiSolverInterface.hpp"
#include "CoinWarmStartBasis.hpp"

#include "BonTMINLP.hpp"
#include "BonTMINLP2TNLP.hpp"
#include "BonTNLP2FPNLP.hpp"
#include "BonTNLPSolver.hpp"
#include "BonCutStrengthener.hpp"
#include "BonRegisteredOptions.hpp"

namespace Bonmin {

  class StrongBranchingSolver;

  /** Solvers for solving nonlinear programs.*/
  enum Solver{
    EIpopt=0 /** <a href="http://projects.coin-or.org/Ipopt">
    Ipopt </a> interior point algorithm.*/,
    EFilterSQP /** <a href="http://www-unix.mcs.anl.gov/~leyffer/solvers.html"> filterSQP </a> Sequential Quadratic Programming algorithm.*/
  };
/**
   This is class provides an Osi interface for a Mixed Integer Linear Program
   expressed as a TMINLP
   (so that we can use it for example as the continuous solver in Cbc).
*/

class OsiTMINLPInterface : public OsiSolverInterface
{
  friend class BonminParam;

public:

  //#############################################################################

  /**Error class to throw exceptions from OsiTMINLPInterface.
   * Inherited from CoinError, we just want to have a different class to be able to catch
   * errors thrown by OsiTMINLPInterface.
  */
class SimpleError : public CoinError
  {
  private:
    SimpleError();

  public:
    ///Alternate constructor using strings
    SimpleError(std::string message,
        std::string methodName,
	std::string f = std::string(),
	int l = -1)
        :
        CoinError(message,methodName,std::string("OsiTMINLPInterface"), f, l)
    {}
  }
  ;

#ifdef __LINE__
#define SimpleError(x, y) SimpleError((x), (y), __FILE__, __LINE__)
#endif

  // Error when problem is not solved
  TNLPSolver::UnsolvedError * newUnsolvedError(int num, Ipopt::SmartPtr<TMINLP2TNLP> problem, std::string name){
    return app_->newUnsolvedError(num, problem, name);
  }
  //#############################################################################


  /** Type of the messages specifically outputed by OsiTMINLPInterface.*/
  enum MessagesTypes{
    SOLUTION_FOUND/**found a feasible solution*/,
    INFEASIBLE_SOLUTION_FOUND/**found an infeasible problem*/,
    UNSOLVED_PROBLEM_FOUND/**found an unsolved problem*/,
    WARNING_RESOLVING /** Warn that a problem is resolved*/,
    WARN_SUCCESS_WS/** Problem not solved with warm start but solved without*/,
    WARN_SUCCESS_RANDOM/** Subproblem not solve with warm start but solved with random point*/,
    WARN_CONTINUING_ON_FAILURE/** a failure occured but is continuing*/,
    SUSPECT_PROBLEM/** Output the number of the problem.*/,
    SUSPECT_PROBLEM2/** Output the number of the problem.*/,
    IPOPT_SUMMARY /** Output summary statistics on Ipopt solution.*/,
    BETTER_SOL /** Found a better solution with random values.*/,
    LOG_HEAD/** Head of "civilized" log.*/,
    LOG_FIRST_LINE/** First line (first solve) of log.*/,
    LOG_LINE/**standard line (retry solving) of log.*/,
    ALTERNATE_OBJECTIVE/** Recomputed integer feasible with alternate objective function*/,
    WARN_RESOLVE_BEFORE_INITIAL_SOLVE /** resolve() has been called but there
                                              was no previous call to initialSolve().
                                         */,
    ERROR_NO_TNLPSOLVER /** Trying to access non-existent TNLPSolver*/,
    WARNING_NON_CONVEX_OA /** Warn that there are equality or ranged constraints and OA may works bad.*/,
    OSITMINLPINTERFACE_DUMMY_END
  };

  //#############################################################################


  /** Messages outputed by an OsiTMINLPInterface. */
class Messages : public CoinMessages
  {
  public:
    /// Constructor
    Messages();
  };


  //#############################################################################


  /**@name Constructors and destructors */
  //@{
  /// Default Constructor
  OsiTMINLPInterface();

  /** Facilitator to initialize interface. */
  void initialize(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
                  Ipopt::SmartPtr<Ipopt::OptionsList> options,
                  Ipopt::SmartPtr<Ipopt::Journalist> journalist_,
                  Ipopt::SmartPtr<TMINLP> tminlp);

  /** Set the model to be solved by interface.*/
  void setModel(Ipopt::SmartPtr<TMINLP> tminlp);
  /** Set the solver to be used by interface.*/
  void setSolver(Ipopt::SmartPtr<TNLPSolver> app);
  /** Sets the TMINLP2TNLP to be used by the interface.*/
  void use(Ipopt::SmartPtr<TMINLP2TNLP> tminlp2tnlp){
     problem_ = tminlp2tnlp;}
  /** Copy constructor
  */
  OsiTMINLPInterface (const OsiTMINLPInterface &);

  /** Virtual copy constructor */
  OsiSolverInterface * clone(bool copyData = true) const;

  /// Assignment operator
  OsiTMINLPInterface & operator=(const OsiTMINLPInterface& rhs);

  /// Destructor
  virtual ~OsiTMINLPInterface ();


  /// Read parameter file
  void readOptionFile(const std::string & fileName);

  /// Retrieve OsiTMINLPApplication option list
  Ipopt::SmartPtr<Ipopt::OptionsList> options();

  //---------------------------------------------------------------------------
  /**@name Solve methods */
  //@{
  /// Solve initial continuous relaxation
  virtual void initialSolve();

  /** Resolve the continuous relaxation after problem modification.
      Have to call initialSolve before calling this
   */
  virtual void resolve();

  /** Resolve the problem with different random starting points
      to try to find a better solution (only makes sense for a non-convex problem.*/
  virtual void resolveForCost(int numretry);

  /** Method to be called when a problem has failed to be solved. Will try
      to resolve it with different settings.
  */
  virtual void resolveForRobustness(int numretry);

  /// Nescessary for compatibility with OsiSolverInterface but does nothing.
  virtual void branchAndBound()
  {
    throw SimpleError("Function not implemented for OsiTMINLPInterface","branchAndBound()");
  }
  //@}



  //---------------------------------------------------------------------------
  ///@name Methods returning info on how the solution process terminated
  //@{
  /// Are there a numerical difficulties?
  virtual bool isAbandoned() const;
  /// Is optimality proven?
  virtual bool isProvenOptimal() const;
  /// Is primal infeasiblity proven?
  virtual bool isProvenPrimalInfeasible() const;
  /// Is dual infeasiblity proven?
  virtual bool isProvenDualInfeasible() const;
  /// Is the given primal objective limit reached?
  virtual bool isPrimalObjectiveLimitReached() const;
  /// Is the given dual objective limit reached?
  virtual bool isDualObjectiveLimitReached() const;
  /// Iteration limit reached?
  virtual bool isIterationLimitReached() const;

  ///Warn solver that branch-and-bound is continuing after a failure
  void continuingOnAFailure()
  {
    hasContinuedAfterNlpFailure_ = true;
  }
  /// Did we continue on a failure
  bool hasContinuedOnAFailure()
  {
    return hasContinuedAfterNlpFailure_;
  }
  /// tell to ignore the failures (don't throw, don't fathom, don't report)
  void ignoreFailures()
  {
    pretendFailIsInfeasible_ = 2;
  }
  /// Force current solution to be infeasible
  void forceInfeasible()
  {
    problem_->set_obj_value(1e200);
  }
  /// Force current solution to be branched on (make it fractionnal with small objective)
  void forceBranchable()
  {
    problem_->set_obj_value(-1e200);
    problem_->force_fractionnal_sol();
  }
  //@}


  //---------------------------------------------------------------------------
  /**@name Parameter set/get methods

     The set methods return true if the parameter was set to the given value,
     false otherwise. There can be various reasons for failure: the given
     parameter is not applicable for the solver (e.g., refactorization
     frequency for the clp algorithm), the parameter is not yet implemented
     for the solver or simply the value of the parameter is out of the range
     the solver accepts. If a parameter setting call returns false check the
     details of your solver.

     The get methods return true if the given parameter is applicable for the
     solver and is implemented. In this case the value of the parameter is
     returned in the second argument. Otherwise they return false.
  */
  //@{
  // Set an integer parameter
  bool setIntParam(OsiIntParam key, int value);
  // Set an double parameter
  bool setDblParam(OsiDblParam key, double value);
  // Set a string parameter
  bool setStrParam(OsiStrParam key, const std::string & value);
  // Get an integer parameter
  bool getIntParam(OsiIntParam key, int& value) const;
  // Get an double parameter
  bool getDblParam(OsiDblParam key, double& value) const;
  // Get a string parameter
  bool getStrParam(OsiStrParam key, std::string& value) const;

  // Get the push values for starting point
  inline double getPushFact() const
  {
    return pushValue_;
  }

  //@}

  /** @name Enums for optionslist parameters */
  //@{
  enum VarSelectStra_Enum {
    MOST_FRACTIONAL=0,
    STRONG_BRANCHING,
    RELIABILITY_BRANCHING,
    CURVATURE_ESTIMATOR,
    QP_STRONG_BRANCHING,
    LP_STRONG_BRANCHING,
    NLP_STRONG_BRANCHING,
    OSI_SIMPLE,
    OSI_STRONG
  };
  //@}

  //---------------------------------------------------------------------------
  /**@name Problem information methods

     These methods call the solver's query routines to return
     information about the problem referred to by the current object.
     Querying a problem that has no data associated with it result in
     zeros for the number of rows and columns, and NULL pointers from
     the methods that return vectors.

     Const pointers returned from any data-query method are valid as
     long as the data is unchanged and the solver is not called.
  */
  //@{
  /**@name Methods related to querying the input data */
  //@{
  /// Get number of columns
  virtual int getNumCols() const;

  /// Get number of rows
  virtual int getNumRows() const;

  ///get name of variables
  const OsiSolverInterface::OsiNameVec& getVarNames() ;
  /// Get pointer to array[getNumCols()] of column lower bounds
  virtual const double * getColLower() const;

  /// Get pointer to array[getNumCols()] of column upper bounds
  virtual const double * getColUpper() const;

  /** Get pointer to array[getNumRows()] of row constraint senses.
      <ul>
      <li>'L': <= constraint
      <li>'E': =  constraint
      <li>'G': >= constraint
      <li>'R': ranged constraint
      <li>'N': free constraint
      </ul>
  */
  virtual const char * getRowSense() const;

  /** Get pointer to array[getNumRows()] of rows right-hand sides
      <ul>
      <li> if rowsense()[i] == 'L' then rhs()[i] == rowupper()[i]
      <li> if rowsense()[i] == 'G' then rhs()[i] == rowlower()[i]
      <li> if rowsense()[i] == 'R' then rhs()[i] == rowupper()[i]
      <li> if rowsense()[i] == 'N' then rhs()[i] == 0.0
      </ul>
  */
  virtual const double * getRightHandSide() const;

  /** Get pointer to array[getNumRows()] of row ranges.
      <ul>
      <li> if rowsense()[i] == 'R' then
      rowrange()[i] == rowupper()[i] - rowlower()[i]
      <li> if rowsense()[i] != 'R' then
      rowrange()[i] is 0.0
      </ul>
  */
  virtual const double * getRowRange() const;

  /// Get pointer to array[getNumRows()] of row lower bounds
  virtual const double * getRowLower() const;

  /// Get pointer to array[getNumRows()] of row upper bounds
  virtual const double * getRowUpper() const;

  /** Get objective function sense (1 for min (default), -1 for max)
   * Always minimizes */
  virtual double getObjSense() const
  {
    return 1;
  }

  /// Return true if column is continuous
  virtual bool isContinuous(int colNumber) const;

  /// Return true if column is binary
  virtual bool isBinary(int columnNumber) const;

  /** Return true if column is integer.
      Note: This function returns true if the the column
      is binary or a general integer.
  */
  virtual bool isInteger(int columnNumber) const;

  /// Return true if column is general integer
  virtual bool isIntegerNonBinary(int columnNumber) const;

  /// Return true if column is binary and not fixed at either bound
  virtual bool isFreeBinary(int columnNumber) const;

  /// Get solver's value for infinity
  virtual double getInfinity() const;

  ///Get priorities on integer variables.
  const int * getPriorities() const
  {
    const TMINLP::BranchingInfo * branch = tminlp_->branchingInfo();
    if(branch)
      return branch->priorities;
    else return NULL;
  }
  ///get prefered branching directions
  const int * getBranchingDirections() const
  {
    const TMINLP::BranchingInfo * branch = tminlp_->branchingInfo();
    if(branch)
      return branch->branchingDirections;
    else return NULL;
  }
  const double * getUpPsCosts() const
  {
    const TMINLP::BranchingInfo * branch = tminlp_->branchingInfo();
    if(branch)
    return branch->upPsCosts;
    else return NULL;
  }
  const double * getDownPsCosts() const
  {
    const TMINLP::BranchingInfo * branch = tminlp_->branchingInfo();
    if(branch)
    return branch->downPsCosts;
    else return NULL;
  }


  //@}

  /**@name Methods related to querying the solution */
  //@{
  /// Get pointer to array[getNumCols()] of primal solution vector
  virtual const double * getColSolution() const;

  /// Get pointer to array[getNumRows()] of dual prices
  virtual const double * getRowPrice() const;

  /// Get a pointer to array[getNumCols()] of reduced costs
  virtual const double * getReducedCost() const;

  /** Get pointer to array[getNumRows()] of row activity levels (constraint
      matrix times the solution vector */
  virtual const double * getRowActivity() const;


  /** Get how many iterations it took to solve the problem (whatever
      "iteration" mean to the solver.
      * \todo Figure out what it could mean for Ipopt.
      */
  virtual int getIterationCount() const;

  /** get total number of calls to solve.*/
  int nCallOptimizeTNLP()
  {
    return nCallOptimizeTNLP_;
  }
  /** get total time taken to solve NLP's. */
  double totalNlpSolveTime()
  {
    return totalNlpSolveTime_;
  }
  /** get total number of iterations */
  int totalIterations()
  {
    return totalIterations_;
  }


  //@}
  //-------------------------------------------------------------------------
  /**@name Methods to modify the objective, bounds, and solution
  */
  //@{

  /** Set a single column lower bound.
      Use -getInfinity() for -infinity. */
  virtual void setColLower( int elementIndex, double elementValue );

  /** Set a single column upper bound.
      Use getInfinity() for infinity. */
  virtual void setColUpper( int elementIndex, double elementValue );

  /** Set the lower bounds for all columns
      array [getNumCols()] is an array of values for the objective.
  */
  virtual void setColLower(const double * array);

  /** Set the upper bounds for all columns
      array [getNumCols()] is an array of values for the objective.
  */
  virtual void setColUpper(const double * array);


  /** Set a single row lower bound.
      Use -getInfinity() for -infinity. */
  virtual void setRowLower( int elementIndex, double elementValue );

  /** Set a single row upper bound.
      Use getInfinity() for infinity. */
  virtual void setRowUpper( int elementIndex, double elementValue );

  /** Set the type of a single row */
  virtual void setRowType(int index, char sense, double rightHandSide,
      double range);


  /** \brief Set the objective function sense (disabled).
   * (1 for min (default), -1 for max)
   \todo Make it work.
   \bug Can not treat maximisation problems. */
  virtual void setObjSense(double s);

  /** Set the primal solution variable values
      Set the values for the starting point.
      \warning getColSolution will never return this vector (unless it is optimal).
  */
  virtual void setColSolution(const double *colsol);

  /** Set dual solution variable values.
      set the values for the starting point.
      \warning getRowPrice will never return this vector (unless it is optimal).
  */
  virtual void setRowPrice(const double * rowprice);

  //@}


  //---------------------------------------------------------------------------
  /**@name WarmStart related methods (those should really do nothing for the moment)*/
  //@{

  /*! \brief Get an empty warm start object

  This routine returns an empty CoinWarmStartBasis object. Its purpose is
  to provide a way to give a client a warm start basis object of the
  appropriate type, which can resized and modified as desired.
  */
  virtual CoinWarmStart *getEmptyWarmStart () const;

  /** Get warmstarting information */
  virtual CoinWarmStart* getWarmStart() const;

  /** Set warmstarting information. Return true/false depending on whether
      the warmstart information was accepted or not. */
  virtual bool setWarmStart(const CoinWarmStart* warmstart);

  void setExposeWarmStart(bool value) {
    exposeWarmStart_ = value;
  }

  void randomStartingPoint();

  //Returns true if a basis is available
  virtual bool basisIsAvailable() const
  {
    // Throw an exception
    throw SimpleError("Needs coding for this interface", "basisIsAvailable");
  }


  //@}

  //-------------------------------------------------------------------------
  /**@name Methods to set variable type */
  //@{
  /** Set the index-th variable to be a continuous variable */
  virtual void setContinuous(int index);
  /** Set the index-th variable to be an integer variable */
  virtual void setInteger(int index);
  //@}

  //Set numIterationSuspect_
  void setNumIterationSuspect(int value)
  {
    numIterationSuspect_ = value;
  }

  /**@name Dummy functions
   * Functions which have to be implemented in an OsiSolverInterface,
   * but which do not do anything (but throwing exceptions) here in the case of a
   * minlp solved using an nlp solver for continuous relaxations */
  //@{

  /** Cbc will understand that no matrix exsits if return -1
  */
  virtual int getNumElements() const
  {
    return -1;
  }


  /** This returns the objective function gradient at the current
   *  point.  It seems to be required for Cbc's pseudo cost
   *  initialization
  */
  virtual const double * getObjCoefficients() const;

  /** We have to keep this but it will return NULL.
   */
  virtual const CoinPackedMatrix * getMatrixByRow() const
  {
      return NULL;
  }


  /** We have to keep this but it will return NULL.
   */
  virtual const CoinPackedMatrix * getMatrixByCol() const
  {
      return NULL;
  }

  /** We have to keep this but it will throw an error.
  */
  virtual void setObjCoeff( int elementIndex, double elementValue )
  {
    throw SimpleError("OsiTMINLPInterface does not implement this function.",
        "setObjCoeff");
  }

  /** We have to keep this but it will throw an error.
  */
  virtual void addCol(const CoinPackedVectorBase& vec,
      const double collb, const double colub,
      const double obj)
  {
    throw SimpleError("OsiTMINLPInterface does not implement this function.",
        "addCol");
  }
  /** We have to keep this but it will throw an error.
  */
  virtual void deleteCols(const int num, const int * colIndices)
  {
    throw SimpleError("OsiTMINLPInterface does not implement this function.",
        "deleteCols");
  }

  /** We have to keep this but it will throw an error.
  */
  virtual void addRow(const CoinPackedVectorBase& vec,
      const double rowlb, const double rowub)
  {
    throw SimpleError("OsiTMINLPInterface does not implement this function.",
        "addRow");
  }
  /** We have to keep this but it will throw an error.
  */
  virtual void addRow(const CoinPackedVectorBase& vec,
      const char rowsen, const double rowrhs,
      const double rowrng)
  {
    throw SimpleError("OsiTMINLPInterface model does not implement this function.",
        "addRow");
  }
  /** We have to keep this but it will throw an error.
  */
  virtual void deleteRows(const int num, const int * rowIndices)
  {
    if(num)
      freeCachedRowRim();
     problem_->removeCuts(num, rowIndices);
  }

  void deleteLastRows(int number){
    if(number)
      freeCachedRowRim();
  problem_->removeLastCuts(number);
  }

  /** We have to keep this but it will throw an error
  */
  virtual void loadProblem(const CoinPackedMatrix& matrix,
      const double* collb, const double* colub,
      const double* obj,
      const double* rowlb, const double* rowub)
  {
    throw SimpleError("OsiTMINLPInterface does not implement this function.",
        "loadProblem");
  }


  /** We have to keep this but it will throw an error.
  */
  virtual void assignProblem(CoinPackedMatrix*& matrix,
      double*& collb, double*& colub, double*& obj,
      double*& rowlb, double*& rowub)
  {
    throw SimpleError("OsiTMINLPInterface does not implement this function.",
        "assignProblem");
  }

  /** We have to keep this but it will throw an error.
  */
  virtual void loadProblem(const CoinPackedMatrix& matrix,
      const double* collb, const double* colub,
      const double* obj,
      const char* rowsen, const double* rowrhs,
      const double* rowrng)
  {
    throw SimpleError("OsiTMINLPInterface does not implement this function.",
        "loadProblem");
  }

  /** We have to keep this but it will throw an error.
  */
  virtual void assignProblem(CoinPackedMatrix*& matrix,
      double*& collb, double*& colub, double*& obj,
      char*& rowsen, double*& rowrhs,
      double*& rowrng)
  {
    throw SimpleError("OsiTMINLPInterface does not implement this function.",
        "assignProblem");
  }


  /** We have to keep this but it will throw an error.
  */
  virtual void loadProblem(const int numcols, const int numrows,
      const int* start, const int* index,
      const double* value,
      const double* collb, const double* colub,
      const double* obj,
      const double* rowlb, const double* rowub)
  {
    throw SimpleError("OsiTMINLPInterface does not implement this function.",
        "loadProblem");
  }

  /** We have to keep this but it will throw an error.
  */
  virtual void loadProblem(const int numcols, const int numrows,
      const int* start, const int* index,
      const double* value,
      const double* collb, const double* colub,
      const double* obj,
      const char* rowsen, const double* rowrhs,
      const double* rowrng)
  {
    throw SimpleError("OsiTMINLPInterface model does not implement this function.",
        "loadProblem");
  }

  /** We have to keep this but it will throw an error.
  */
  virtual int readMps(const char *filename,
      const char *extension = "mps")
  {
    throw SimpleError("OsiTMINLPInterface does not implement this function.",
        "readMps");
  }


  /** We have to keep this but it will throw an error.
  */
  virtual void writeMps(const char *filename,
      const char *extension = "mps",
      double objSense=0.0) const
  {
    throw SimpleError("OsiTMINLPInterface does not implement this function.",
        "writeMps");
  }

  /** Throws an error */
  virtual std::vector<double*> getDualRays(int maxNumRays) const
  {
    throw SimpleError("OsiTMINLPInterface does not implement this function.",
        "getDualRays");
  }

  /** Throws an error */
  virtual std::vector<double*> getPrimalRays(int maxNumRays) const
  {
    throw SimpleError("OsiTMINLPInterface does not implement this function.",
        "getPrimalRays");
  }

  //@}


  
  //---------------------------------------------------------------------------



  /**@name Control of Ipopt output
   */
  //@{
  void turnOffSolverOutput(){
  app_->turnOffOutput();}
  void turnOnSolverOutput(){
  app_->turnOnOutput();}
  //@}

  /**@name Sets and Getss */
  //@{
  /// Get objective function value (can't use default)
  virtual double getObjValue() const;

  //@}

  /** get pointer to the TMINLP2TNLP adapter */
  const TMINLP2TNLP * problem() const
  {
    return GetRawPtr(problem_);
  }

  TMINLP2TNLP * problem()
  {
    return GetRawPtr(problem_);
  }

  const TMINLP * model() const
  {
    return GetRawPtr(tminlp_);
  }
  
  Bonmin::TMINLP * model()
  {
    return GetRawPtr(tminlp_);
  }
  
  const Bonmin::TNLPSolver * solver() const
  {
    return GetRawPtr(app_);
  } 
 
  TNLPSolver * solver()
  {
    return GetRawPtr(app_);
  } 
  /** \name Methods to build outer approximations */
  //@{
  /** \name Methods to build outer approximations */
  //@{
  /** \brief Extract a linear relaxation of the MINLP.
   * Use user-provided point to build first-order outer-approximation constraints at the optimum.
   * And put it in an OsiSolverInterface.
   */
  virtual void extractLinearRelaxation(OsiSolverInterface &si, const double *x, 
                                       bool getObj = 1);

  /** \brief Extract a linear relaxation of the MINLP.
   * Solve the continuous relaxation and takes first-order outer-approximation constraints at the optimum.
   * The put everything in an OsiSolverInterface.
   */
  virtual void extractLinearRelaxation(OsiSolverInterface &si, bool getObj = 1,
                                       bool solveNlp = 1){
     if(solveNlp)
       initialSolve();
     extractLinearRelaxation(si, getColSolution(), getObj); 
     if(solveNlp){
        app_->enableWarmStart();
        setColSolution(problem()->x_sol());
        setRowPrice(problem()->duals_sol());
     }
   }

  /** Get the outer approximation constraints at the current optimal point.
      If x2 is different from NULL only add cuts violated by x2.
   (Only get outer-approximations of nonlinear constraints of the problem.)*/
  void getOuterApproximation(OsiCuts &cs, bool getObj, const double * x2, bool global)
{
  getOuterApproximation(cs, getColSolution(), getObj, x2, global);
}

  /** Get the outer approximation constraints at provided point.
      If x2 is different from NULL only add cuts violated by x2.
   (Only get outer-approximations of nonlinear constraints of the problem.)*/
  void getOuterApproximation(OsiCuts &cs, const double * x, bool getObj, const double * x2, bool global){
    getOuterApproximation(cs, x, getObj, x2, 0., global);}

  /** Get the outer approximation constraints at provided point.
      If x2 is different from NULL only add cuts violated by x2 by more than delta.
   (Only get outer-approximations of nonlinear constraints of the problem.)*/
  virtual void getOuterApproximation(OsiCuts &cs, const double * x, bool getObj, const double * x2,
                                     double theta, bool global);

 /** Get the outer approximation at provided point for given constraint. */
  virtual void getConstraintOuterApproximation(OsiCuts & cs, int constraintNumber,
                                               const double * x, 
                                               const double * x2, bool global);

 /** Get the outer approximation at current optimal point for given constraint. */
  void getConstraintOuterApproximation(OsiCuts & cs, int constraintNumber,
                                       const double * x2, bool global){
     getConstraintOuterApproximation(cs, constraintNumber, getColSolution(),x2,global);
  }

  /** Get the Benders cut at provided point with provided multipliers.*/
  void getBendersCut(OsiCuts &cs, const double * x, const double *lambda, bool getObj = 1);

  /** solve the problem of finding the closest point to x_bar in the subspace of coordinates given by ind
   * (i.e., \f$ min \sum\limits_{i=1}^n (x_{ind[i]} -\overline{x}_i)^2 \f$ ,
   * and get the corresponding outer-approximation constraints.
      (Only get outer-approximations of nonlinear constraints of the problem.)
   * \return Distance between feasibility set and x
   * \param n number of element in arrays x and ind
   * \param ind indices of the coordinate*/
  double getFeasibilityOuterApproximation(int n, const double * x_bar,const int *ind, OsiCuts &cs, bool addOnlyViolated, bool global);
  //@}

  /** \name output for OA cut generation
       \todo All OA code here should be moved to a separate class sometime.*/
  //@{
  /** OA Messages types.*/
  enum OaMessagesTypes {
    CUT_NOT_VIOLATED_ENOUGH = 0/** Says that one cut has been generarted, where from, which is the violation.*/,
    VIOLATED_OA_CUT_GENERATED/** Cut is not violated enough, give violation.*/,
    OA_CUT_GENERATED/** Print the cut which has been generated.*/,
    OA_MESSAGES_DUMMY_END/** Dummy end.*/};
  /** Class to store OA Messages.*/
  class OaMessages :public CoinMessages{
     public:
     /** Default constructor.*/
     OaMessages();
  };
  /** Like a CoinMessageHandler but can print a cut also.*/
  class OaMessageHandler : public CoinMessageHandler{
    public:
    /** Default constructor.*/
    OaMessageHandler():CoinMessageHandler(){
    }
    /** Constructor to put to file pointer (fp won't be closed).*/
    OaMessageHandler(FILE * fp):CoinMessageHandler(fp){
    }
    /** Destructor.*/
    virtual ~OaMessageHandler(){
    }
    /** Copy constructor.*/
    OaMessageHandler(const OaMessageHandler &other):
    CoinMessageHandler(other){}
    /** Constructor from a regular CoinMessageHandler.*/
    OaMessageHandler(const CoinMessageHandler &other):
    CoinMessageHandler(other){}
    /** Assignment operator.*/
    OaMessageHandler & operator=(const OaMessageHandler &rhs){
       CoinMessageHandler::operator=(rhs);
       return *this;}
    /** Virtual copy */
    virtual CoinMessageHandler* clone() const{
      return new OaMessageHandler(*this);}
    /** print an OsiRowCut.*/
    void print(OsiRowCut &row);
  };
  void setOaMessageHandler(const CoinMessageHandler &handler){
    delete oaHandler_;
    oaHandler_ = new OaMessageHandler(handler);
  }
  //@}

    //-----------------------------------------------------------------------
    /** Apply a collection of cuts.
    */
    virtual ApplyCutsReturnCode applyCuts(const OsiCuts & cs,
					  double effectivenessLb = 0.0){
      problem_->addCuts(cs);
      ApplyCutsReturnCode rc;
      return rc;}

   /** Add a collection of linear cuts to problem formulation.*/
  virtual void applyRowCuts(int numberCuts, const OsiRowCut * cuts);


  /** Add a collection of linear cuts to the problem formulation */
  virtual void applyRowCuts(int numberCuts, const OsiRowCut ** cuts)
  {
    if(numberCuts)
      freeCachedRowRim();
    problem_->addCuts(numberCuts, cuts);
  }

 /** Get infinity norm of constraint violation for x. Put into
     obj the objective value of x.*/
 double getConstraintsViolation(const double * x, double & obj);

  /** Get infinity norm of constraint violation for x and error in objective
      value where obj is the estimated objective value of x.*/
  double getNonLinearitiesViolation(const double *x, const double obj);

//---------------------------------------------------------------------------

  void extractInterfaceParams();


  /** To set some application specific defaults. */
  virtual void setAppDefaultOptions(Ipopt::SmartPtr<Ipopt::OptionsList> Options);

  /** Register all possible options to Bonmin */
  static void registerOptions (Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions);
  
  Ipopt::SmartPtr<Bonmin::RegisteredOptions> regOptions(){
    if(IsValid(app_))
      return app_->roptions();
    else
      return NULL;
  }

  /** @name Methods related to strong branching */
  //@{
  /// Set the strong branching solver
  void SetStrongBrachingSolver(Ipopt::SmartPtr<StrongBranchingSolver> strong_branching_solver);
  /// Create a hot start snapshot of the optimization process.  In our
  /// case, we initialize the StrongBrachingSolver.
  virtual void markHotStart();
  /// Optimize starting from the hot start snapshot. In our case, we
  /// call the StrongBranchingSolver to give us an approximate
  /// solution for the current state of the bounds
  virtual void solveFromHotStart();
  /// Delete the hot start snapshot. In our case we deactivate the
  /// StrongBrachingSolver.
  virtual void unmarkHotStart();
  //@}

protected:
  
  //@}

  enum RandomGenerationType{
    uniform =0, perturb=1, perturb_suffix=2};
  /// Initialize data structures for storing the jacobian
  int initializeJacobianArrays();

  ///@name Virtual callbacks for application specific stuff
  //@{
  virtual std::string  appName()
  {
    return "bonmin";
  }
  //@}
  ///@name Protected methods
  //@{

  /** Call Ipopt to solve or resolve the problem and check for errors.*/
  void solveAndCheckErrors(bool doResolve, bool throwOnFailure,
      const char * whereFrom);


  /** Add a linear cut to the problem formulation.
  */
  virtual void applyRowCut( const OsiRowCut & rc )
  {
    const OsiRowCut * cut = &rc;
    problem_->addCuts(1, &cut);
  }
  /** We have to keep this but it will throw an error.
  */
  virtual void applyColCut( const OsiColCut & cc )
  {
    throw SimpleError("Ipopt model does not implement this function.",
        "applyColCut");
  }

//  /** Read the name of the variables in an ampl .col file. */
//  void readVarNames() const;

  //@}

  /**@name Model and solver */
  //@{
  /** TMINLP model.*/
  Ipopt::SmartPtr<TMINLP> tminlp_;
  /** Adapter for a MINLP to a NLP */
  Ipopt::SmartPtr<TMINLP2TNLP> problem_;
  /** Solver for a TMINLP. */
  Ipopt::SmartPtr<TNLPSolver> app_;
  //@}

  /**@name Cached information on the problem */
  //@{
  /** Free cached data relative to variables */
  void freeCachedColRim();
  /** Free cached data relative to constraints */
  void freeCachedRowRim();
  /** Free all cached data*/
  void freeCachedData();
  /** Extract rowsense_ vector rhs_ vector and rowrange_ vector from the lower and upper bounds
   *  on the constraints */
  void extractSenseRhsAndRange() const;
  /// Pointer to dense vector of row sense indicators
  mutable char    *rowsense_;

  /// Pointer to dense vector of row right-hand side values
  mutable double  *rhs_;

  /// Pointer to dense vector of slack upper bounds for range constraints (undefined for non-range rows)
  mutable double  *rowrange_;
  /** Pointer to dense vector of reduced costs
      \warning Always 0. with Ipopt*/
  mutable double  *reducedCosts_;
  /** DualObjectiveLimit is used to store the cutoff in Cbc*/
  double OsiDualObjectiveLimit_;
  /** does the file variable names exists (will check automatically).*/
  mutable bool hasVarNamesFile_;
  //@}
  /// number of time NLP has been solved
  int nCallOptimizeTNLP_;
  /// Total solution time of NLP
  double totalNlpSolveTime_;
  /// toatal number of iterations
  int totalIterations_;
  /// max radius for random point
  double maxRandomRadius_;
  /// Method to pick a random starting point.
  int randomGenerationType_;
  /// Maximum perturbation value
  double max_perturbation_;
  /// Ipopt value for pushing initial point inside the bounds
  double pushValue_;
  /// Number of times problem will be resolved in initialSolve (root node)
  int numRetryInitial_;
  /// Number of times problem will be resolved in resolve
  int numRetryResolve_;
  /// Number of times infeasible problem will be resolved.
  int numRetryInfeasibles_;
  /// Number of times problem will be resolved in case of a failure
  int numRetryUnsolved_;
  /** Messages specific to an OsiTMINLPInterface. */
  Messages messages_;
  /** If not 0 when a problem is not solved (failed to be solved)
      will pretend that it is infeasible. If == 1 will care
      (i.e. record the fact issue messages to user), if ==2 don't care (somebody else will) */
  int pretendFailIsInfeasible_;
  /** did we ever continue optimization ignoring a failure. */
  bool hasContinuedAfterNlpFailure_;
  /** number iterations above which a problem is considered suspect (-1 is considered \f$+ \infty \f$).
  	If in a call to solve a problem takes more than that number of iterations it will be outputed to files.*/
  int numIterationSuspect_ ;
  /** Has problem been optimized since last change (include setColSolution).
     If yes getColSolution will return Ipopt point, otherwise will return
     initial point.*/
  bool hasBeenOptimized_;
  /** A fake objective function (all variables to 1) to please Cbc
      pseudo costs initialization.  AW: I changed this, it will now be
      the objective gradient at current point. */
  mutable double * obj_;
  /** flag to say wether options have been printed or not.*/
  static bool hasPrintedOptions;

  /** Adapter for TNLP to a feasibility problem */
  Ipopt::SmartPtr<TNLP2FPNLP> feasibilityProblem_;


  /** \name Arrays to store Jacobian matrix */
  //@{
  /** Row indices.*/
  int * jRow_;
  /** Column indices.*/
  int * jCol_;
  /** Values */
  double * jValues_;
  /** Number of elements.*/
  int nnz_jac;
  //@}

  ///Store the types of the constraints (linear and nonlinear).
  Ipopt::TNLP::LinearityType * constTypes_;
  /* Numerotation of linear/nonlinear constraints
   * Perform independent numerotation of linear (resp. nonlinear constraints)
   * so that constraints of each type are numeroted consecutively */
 // int * constTypesNum_;
  /** Number of linear constraints */
  int nLinear_;
  /** Number of nonlinear constraint
   */
  int nNonLinear_;
  /** Value for small non-zero element which we will try to remove cleanly in OA cuts.*/
  double tiny_;
  /** Value for small non-zero element which we will take the risk to ignore in OA cuts.*/
  double veryTiny_;
  /** Value for infinity. */
  double infty_;
  /** status of last optimization. */
  TNLPSolver::ReturnStatus optimizationStatus_;
  /** Flag indicating if the warm start methods actually do something.*/
  bool exposeWarmStart_;
  /** Is it the first solve (for random starting point at root options).*/
  bool firstSolve_;
  /** Object for strengthening cuts */
  SmartPtr<CutStrengthener> cutStrengthener_;

  /** \name output for OA cut generation
       \todo All OA code here should be moved to a separate class sometime.*/
  //@{
  /** OA Messages.*/
  OaMessages oaMessages_;
  /** OA Message handler. */
  OaMessageHandler * oaHandler_;
  //@}
protected:
  /** Facilitator to create an application. */
  void createApplication(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
                         Ipopt::SmartPtr<Ipopt::OptionsList> options,
                         Ipopt::SmartPtr<Ipopt::Journalist> journalist);
  ///Constructor without model only for derived classes
  OsiTMINLPInterface(Ipopt::SmartPtr<TNLPSolver> app);

private:
  /** solver to be used for all strong branching solves */
  SmartPtr<StrongBranchingSolver> strong_branching_solver_;
  /** status of last optimization before hot start was marked. */
  TNLPSolver::ReturnStatus optimizationStatusBeforeHotStart_;

  // DELETEME
  SmartPtr<StrongBranchingSolver> strong_branching_solver_compare_;
};
}
#endif

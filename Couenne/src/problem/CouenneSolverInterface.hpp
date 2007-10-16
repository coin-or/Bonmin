/*
 * Name:    CouenneSolverInterface.hpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: OsiSolverInterface with a pointer to a CouenneCutGenerator object
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNESOLVERINTERFACE_HPP
#define COUENNESOLVERINTERFACE_HPP

#include "CoinFinite.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CouenneCutGenerator.hpp"


/// Solver interface class with a pointer to a Couenne cut
/// generator. Its main purposes are:
///
/// 1) to apply bound tightening before re-solving
/// 2) to replace OsiSolverInterface::isInteger () with problem_ -> [expression] -> isInteger ()
/// 3) to use NLP solution at branching
 
class CouenneSolverInterface: public OsiClpSolverInterface {

public:

  /// Constructor
  CouenneSolverInterface (CouenneCutGenerator *cg = NULL);

  /// Copy constructor
  CouenneSolverInterface (const CouenneSolverInterface &src);

  /// Destructor
  ~CouenneSolverInterface () {}

  /// Clone
  virtual OsiSolverInterface * clone (bool copyData = true) const
    {return new CouenneSolverInterface (*this);}

  /// we need to overwrite this since we might have internal knowledge
  virtual bool isProvenPrimalInfeasible() const;

  /// we need to overwrite this since we might have internal knowledge
  virtual bool isProvenOptimal() const;

  /// Return cut generator pointer
  CouenneCutGenerator *CutGen ()
    {return cutgen_;}

  /// Set cut generator pointer after setup, to avoid changes in the
  /// pointer due to cut generator cloning (it happens twice in the
  /// algorithm)
  void setCutGenPtr (CouenneCutGenerator *cg)
    {cutgen_ = cg;}

  /// Solve initial LP relaxation 
  virtual void initialSolve (); 

  /// Resolve an LP relaxation after problem modification
  virtual void resolve ();

  /// Resolve an LP without applying bound tightening beforehand
  virtual void resolve_nobt ()
    {OsiClpSolverInterface::resolve ();}

  /** @name Methods for strong branching.
   */
  //@{
  /// Create a hot start snapshot of the optimization process.
  virtual void markHotStart();
  /// Optimize starting from the hot start snapshot.
  virtual void solveFromHotStart();
  /// Delete the hot start snapshot.
  virtual void unmarkHotStart();
  //@}

private:
  /// The pointer to the Couenne cut generator. Gives us a lot of
  /// information, for instance the nlp solver pointer, and the chance
  /// to do bound tightening before resolve ().
  CouenneCutGenerator *cutgen_;

  /// Flag indicating that infeasibility was detected during solveFromHotStart
  bool knowInfeasible_;

  /// Flag indicating that optimality was detected during solveFromHotStart
  bool knowOptimal_;
};

#endif

// (C) Copyright International Business Machines (IBM) 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, International Business Machines
//
// Date :  12/07/2006


// Code separated from BonOaDecBase to try to clarify OAs
#ifndef BonSubMipSolver_HPP
#define BonSubMipSolver_HPP
#include "IpSmartPtr.hpp"
/* forward declarations.*/
class OsiSolverInterface;
class OsiClpSolverInterface;
class OsiCpxSolverInterface;
class CbcStrategy;
class CbcStrategyDefault;
class CbcModel;

namespace Bonmin {
    class RegisteredOptions;
    /** A very simple class to provide a common interface for solving MIPs with Cplex and Cbc.*/
    class SubMipSolver
    {
    public:
      /** Constructor */
      SubMipSolver(OsiSolverInterface * lp = 0,
          const CbcStrategy * strategy = 0);

      ~SubMipSolver();

      /** Assign lp solver. */
      void setLpSolver(OsiSolverInterface * lp);

      /** Assign a strategy. */
      void setStrategy(CbcStrategyDefault * strategy);

      /** get the solution found in last local search (return NULL if no solution). */
      const double * getLastSolution()
      {
        return integerSolution_;
      }

      double getLowerBound()
      {
        return lowBound_;
      }
      /** update cutoff and perform a local search to a good solution. */
      void find_good_sol(double cutoff,
          int loglevel,
          double maxTime);

      /** update cutoff and optimize MIP. */
      void optimize(double cutoff,
          int loglevel,
          double maxTime);

      /** Returns lower bound. */
      inline double lowBound()
      {
        return lowBound_;
      }

      /** returns optimality status. */
      inline bool optimal()
      {
        return optimal_;
      }

      /** Returns number of nodes in last solve.*/
      inline int nodeCount()
      {
        return nodeCount_;
      }

      /** Returns number of simplex iterations in last solve.*/
      inline int iterationCount()
      {
        return iterationCount_;
      }


      /** Register options for that Oa based cut generation method. */
      OsiSolverInterface * solver(){
         return lp_;
      }
     static void registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions);
    private:
      /** lp (potentially mip solver). */
      OsiSolverInterface * lp_;
      /** If lp solver is clp (then have to use Cbc).*/
      OsiClpSolverInterface *clp_;
      /** If lp solver is cpx this points to it. */
      OsiCpxSolverInterface * cpx_;
      /** If cbc is used pointer to CbcModel. */
      CbcModel * cbc_;
      /** lower bound obtained */
      double lowBound_;
      /** Is optimality proven? */
      bool optimal_;
      /** Has an integer solution? then it is here*/
      double * integerSolution_;
      /** Strategy for solving sub mips with cbc. */
      CbcStrategyDefault * strategy_;
      /** number of nodes in last mip solved.*/
      int nodeCount_;
      /** number of simplex iteration in last mip solved.*/
      int iterationCount_;
    };

}

#endif


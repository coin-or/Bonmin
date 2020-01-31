// (C) Copyright International Business Machines (IBM) 2006
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// P. Bonami, International Business Machines
//
// Date :  12/07/2006


// Code separated from BonOaDecBase to try to clarify OAs
// A number of utility to manipulate models used in OA and FP
#ifndef BonSolverHelp_H
#define BonSolverHelp_H

class OsiSolverInterface;
class OsiBranchingInformation;
class OsiObject;
class OsiCuts;

namespace Bonmin {
   /** Check for integer feasibility of a solution return true if it is feasible.*/
   bool integerFeasible(OsiSolverInterface & si, const OsiBranchingInformation & info,
                        double integer_tolerance,
                        OsiObject ** objects = 0, int nObjects = -1);

   /** Fix integer variables in si to their values in colsol.
       \remark colsol is assumed to be integer on the integer constrained variables.
    */
   void fixIntegers(OsiSolverInterface & si, const OsiBranchingInformation & info,
                    double integer_tolerance,
                    OsiObject ** objects = 0, int nObjects = -1);
   /** Relax integer variables in si.
    */
   void relaxIntegers(OsiSolverInterface & si, const OsiBranchingInformation & info,
                    double integer_tolerance,
                    OsiObject ** objects = 0, int nObjects = -1);
   /** Check if two solutions are the same on integer variables. */
   bool isDifferentOnIntegers(OsiSolverInterface &si,
                              OsiObject ** objects, int nObjects,
                              double integer_tolerance,
                              const double * colsol, const double * other);

   /** Install cuts in solver. */
   void installCuts(OsiSolverInterface &si,
                    const OsiCuts& cs, int numberCuts);

}
#endif


/*
 * Name:    CouenneChooseVariable.cpp
 * Authors: Pierre Bonami, IBM Corp.
 *          Pietro Belotti, Carnegie Mellon University
 * Purpose: Branching object for choosing branching auxiliary variable
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <CouenneChooseVariable.hpp>

/// Default Constructor 
CouenneChooseVariable::CouenneChooseVariable (): 
  problem_ (NULL) {}

/// Constructor from solver (so we can set up arrays etc)
CouenneChooseVariable::CouenneChooseVariable (const OsiSolverInterface *,
					      CouenneProblem *p):
  problem_ (p) {}

/// Copy constructor 
CouenneChooseVariable::CouenneChooseVariable (const CouenneChooseVariable &source):
  OsiChooseVariable (source),
  problem_ (source.Problem ()) {}

/// Assignment operator 
CouenneChooseVariable & CouenneChooseVariable::operator= (const CouenneChooseVariable& rhs)
{}

/// Clone
OsiChooseVariable *CouenneChooseVariable::clone() const
  {return new CouenneChooseVariable (*this);}

/// Destructor 
CouenneChooseVariable::~CouenneChooseVariable () {

}

/** Sets up strong list and clears all if initialize is true.
    Returns number of infeasibilities. 
    If returns -1 then has worked out node is infeasible!
*/

int CouenneChooseVariable::setupList (OsiBranchingInformation *info, bool initialize) {

  /*  std::cout<<"Branching on :"<<std::endl;
  for(int kk = 0 ; kk < info -> numberColumns_ ; kk++)
    {
            std::cout<<"x["<<kk<<"] = "<<info->solution_[kk]<<"\t";  
    }
  std::cout<<std::endl;
  */
  problem_ -> update ((CouNumber *) (info -> solution_), 
		      (CouNumber *) (info -> lower_), 
		      (CouNumber *) (info -> upper_));

  return OsiChooseVariable::setupList (info, initialize);
}


/** Choose a variable
    Returns:
    -1 Node is infeasible
    0  Normal termination - we have a candidate
    1  All looks satisfied - no candidate
    2  We can change the bound on a variable - but we also have a strong branching candidate
    3  We can change the bound on a variable - but we have a non-strong branching candidate
    4  We can change the bound on a variable - no other candidates
    We can pick up branch from bestObjectIndex() and bestWhichWay()
    We can pick up a forced branch (can change bound) from firstForcedObjectIndex() 
    and firstForcedWhichWay()
    If we have a solution then we can pick up from goodObjectiveValue() and goodSolution()
    If fixVariables is true then 2,3,4 are all really same as problem changed
*/

int CouenneChooseVariable::chooseVariable (OsiSolverInterface * solver, 
					   OsiBranchingInformation *info, 
					   bool fixVariables) {
  if (numberUnsatisfied_) {
    bestObjectIndex_ = list_ [0];
    bestWhichWay_ = solver -> object (bestObjectIndex_) -> whichWay();
    firstForcedObjectIndex_ = -1;
    firstForcedWhichWay_    = -1;
    return 0;
  } else return 1;
}

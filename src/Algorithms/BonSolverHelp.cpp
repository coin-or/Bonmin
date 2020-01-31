// (C) Copyright International Business Machines (IBM) 2006
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// P. Bonami, International Business Machines
//
// Date :  12/07/2006


// Code separated from BonOaDecBase to try to clarify OAs

#include <sstream>
#include "BonSolverHelp.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiBranchingObject.hpp"
#include "OsiCuts.hpp"
#include "CoinWarmStartBasis.hpp"
namespace Bonmin {
/// Check for integer feasibility of a solution return 1 if it is
bool integerFeasible(OsiSolverInterface & si, const OsiBranchingInformation & info, 
                     double integer_tolerance,
                     OsiObject ** objects, int nObjects) 
{
  if (objects) {
    int dummy;
    for (int i = 0 ; i < nObjects ; i++) {
      double infeasibility = objects[i]->infeasibility(&info, dummy);
      if (infeasibility > 1000*integer_tolerance) return false;
    }
  }
  else {
    const double * sol = info.solution_;
    int numcols = si.getNumCols();
    for (int i = 0 ; i < numcols ; i++) {
      if (si.isInteger(i)) {
        if (fabs(sol[i] - floor(sol[i] + 0.5)) >
            integer_tolerance) {
          return false;
        }
      }
    }
  }
  return true;
}

/** Fix integer variables in si to their values in colsol.
*/
void fixIntegers(OsiSolverInterface & si, 
                 const OsiBranchingInformation & info,
                 double integer_tolerance,
                 OsiObject ** objects, int nObjects)
{
  if (objects) {
    for (int i = 0 ; i < nObjects ; i++) {
      objects[i]->feasibleRegion(&si, &info);
    }
  }
  else {
    const double * colsol = info.solution_;
    for (int i = 0; i < info.numberColumns_; i++) {
      if (si.isInteger(i)) {
        double  value =  colsol[i];
#ifndef NDEBUG
        if (fabs(value - floor(value+0.5)) > integer_tolerance) {
          std::stringstream stream;
          stream<<"Error not integer valued solution"<<std::endl;
          stream<<"---------------- x["<<i<<"] = "<<value<<std::endl;
          throw CoinError(stream.str(),"fixIntegers","OaDecompositionBase::solverManip");
        }
#endif
        value = floor(value+0.5);
        if (fabs(value) > 1e10) {
          std::stringstream stream;
          stream<<"Can not fix variable in nlp because it has too big a value ("<<value
          <<") at optimium of LP relaxation. You should try running the problem with B-BB"<<std::endl;
          throw CoinError(stream.str(),
              "fixIntegers","OaDecompositionBase::solverManip") ;
        }
#ifdef OA_DEBUG
        //         printf("xx %d at %g (bounds %g, %g)",i,value,nlp_->getColLower()[i],
        //                nlp_->getColUpper()[i]);
        std::cout<<(int)value;
#endif
        si.setColLower(i,value);
        si.setColUpper(i,value);
      }
    }
#ifdef OA_DEBUG
    std::cout<<std::endl;
#endif
  }
}

/** Slightly relax integer variables in si.
*/
void relaxIntegers(OsiSolverInterface & si, 
                 const OsiBranchingInformation & info,
                 double integer_tolerance,
                 OsiObject ** objects, int nObjects)
{
  if (objects) {
    for (int i = 0 ; i < nObjects ; i++) {
      OsiSimpleInteger * obj = dynamic_cast<OsiSimpleInteger *>(objects[i]);
      int colNumber = obj->columnNumber();
      si.setColLower(colNumber, si.getColLower()[colNumber] - integer_tolerance); 
      si.setColUpper(colNumber, si.getColUpper()[colNumber] + integer_tolerance); 
    }
  }
  else {
    for (int i = 0; i < info.numberColumns_; i++) {
      if (si.isInteger(i)) {
        const int &colNumber = i;
        si.setColLower(colNumber, si.getColLower()[colNumber] - integer_tolerance); 
        si.setColUpper(colNumber, si.getColUpper()[colNumber] + integer_tolerance); 
      }
    }
  }
}

bool refixIntegers(OsiSolverInterface & si, 
                 const OsiBranchingInformation & info,
                 double integer_tolerance,
                 OsiObject ** objects, int nObjects)
{
  if(!si.isProvenOptimal()) return false;
  if (objects) {
    for (int i = 0 ; i < nObjects ; i++) {
      OsiSimpleInteger * obj = dynamic_cast<OsiSimpleInteger *>(objects[i]);
      int colNumber = obj->columnNumber();
      si.setColLower(colNumber, si.getColLower()[colNumber] - integer_tolerance); 
      si.setColUpper(colNumber, si.getColUpper()[colNumber] + integer_tolerance); 
    }
  }
  else {
    for (int i = 0; i < info.numberColumns_; i++) {
      if (si.isInteger(i)) {
        const int &colNumber = i;
        si.setColLower(colNumber, si.getColLower()[colNumber] - integer_tolerance); 
        si.setColUpper(colNumber, si.getColUpper()[colNumber] + integer_tolerance); 
      }
    }
  }
  return true;
}

/** Install cuts in solver. */
void installCuts(OsiSolverInterface &si,
                 const OsiCuts& cs, int numberCuts){
  int numberCutsBefore = cs.sizeRowCuts() - numberCuts;

  CoinWarmStartBasis * basis
  = dynamic_cast<CoinWarmStartBasis*>(si.getWarmStart()) ;
  assert(basis != NULL); // make sure not volume
  int numrows = si.getNumRows();
  basis->resize(numrows + numberCuts,si.getNumCols());
  for (int i = 0 ; i < numberCuts ; i++) {
    basis->setArtifStatus(numrows + i,
        CoinWarmStartBasis::basic) ;
  }

  const OsiRowCut ** addCuts = new const OsiRowCut * [numberCuts] ;
  for (int i = 0 ; i < numberCuts ; i++) {
    addCuts[i] = &cs.rowCut(i + numberCutsBefore) ;
  }
  si.applyRowCuts(numberCuts,addCuts) ;
  delete [] addCuts ;
  if (si.setWarmStart(basis) == false) {
    delete basis;
    throw CoinError("Fail setWarmStart() after cut installation.",
                    "generateCuts","OACutGenerator2") ;
  }
  delete basis;
}


/** Check if two solutions are the same on integer variables. */
bool isDifferentOnIntegers(OsiSolverInterface &si,
                           OsiObject ** objects, int nObjects,
                           double integer_tolerance,
                           const double * colsol, const double *otherSol)
{
  if (objects) {
    for (int i = 0 ; i < nObjects ; i++) {
      OsiObject * obj = objects[i];
      int colnum = obj->columnNumber();
      if (colnum >= 0) {//Variable branching object
        if (fabs(otherSol[colnum] - colsol[colnum]) > integer_tolerance) {
          return true;
        }
      }
      else {//It is a sos branching object
        OsiSOS * sos = dynamic_cast<OsiSOS *>(obj);
        assert(sos);
        const int * members = sos->members();
        int end = sos->numberMembers();
        for (int k = 0 ; k < end ; k++) {
          if (fabs(otherSol[members[k]] - colsol[members[k]]) > integer_tolerance) {
            return true;
          }
        }
      }
    }
  }
  else {
    int numcols = si.getNumCols();
    for (int i = 0; i < numcols ; i++) {
      if (si.isInteger(i) && fabs(otherSol[i] - colsol[i])>integer_tolerance)
        return true;
    }
  }
  return false;
}

}


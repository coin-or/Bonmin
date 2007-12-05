// (C) Copyright International Business Machines Corporation (IBM) 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pietro Belotti, Carnegie Mellon University
// Pierre Bonami, International Business Machines Corporation
//
// Date : 12/19/2006


#include "BonCouenneInterface.hpp"
#include "CoinHelperFunctions.hpp"
#include <CouenneProblem.hpp>

namespace Bonmin {

/** Default constructor. */
CouenneInterface::CouenneInterface():
  AmplInterface()
{}

/** Copy constructor. */
CouenneInterface::CouenneInterface(const CouenneInterface &other):
  AmplInterface(other)
  {
  }

/** virutal copy constructor. */
CouenneInterface * CouenneInterface::clone(bool CopyData){
  return new CouenneInterface(*this);
}

/** Destructor. */
CouenneInterface::~CouenneInterface(){
}



void 
CouenneInterface::readAmplNlFile(char **& argv, Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
                                 Ipopt::SmartPtr<Ipopt::OptionsList> options,
                                 Ipopt::SmartPtr<Ipopt::Journalist> journalist){
  AmplInterface::readAmplNlFile(argv, roptions, options, journalist);
}

/** \name Overloaded methods to build outer approximations */
  //@{
  /** \brief Extract a linear relaxation of the MINLP.
   *
   * Solve the continuous relaxation and takes first-order
   * outer-approximation constraints at the optimum.  Then put
   * everything in an OsiSolverInterface.
   *
   * The OsiSolverInterface si is empty and has to be populated with
   * the initial linear relaxation.
   */

void 
CouenneInterface::extractLinearRelaxation 
(OsiSolverInterface &si, CouenneCutGenerator & couenneCg, bool getObj, bool solveNlp) {

   if (solveNlp) {

     initialSolve ();

     // TODO: solve nlp, round as in BonNlpHeuristic, and re-solve
     // TODO: better test is checkNLP

     if (getNumIntegers () == 0) { // only if no integer variables

       const CouNumber 
	 obj      = getObjValue ();
	 //	 *solution = getColSolution ();

       // set cutoff to take advantage of bound tightening
       couenneCg. Problem () -> setCutOff (obj);
       /*
       OsiAuxInfo * auxInfo = si. getAuxiliaryInfo ();
       BabInfo * babInfo = dynamic_cast <BabInfo *> (auxInfo);

       if (babInfo) {
	 babInfo -> setNlpSolution (solution, getNumCols (), obj);
	 babInfo -> setHasNlpSolution (true);
       }
       */
     }
   }

   int numcols     = getNumCols (),         // # original               variables
       numcolsconv = couenneCg.getnvars (); // # original + # auxiliary variables

   const double
     *lb = getColLower (),
     *ub = getColUpper ();

   // add original variables to the new problem
   for (register int i=0; i<numcols; i++)
     si.addCol (0, NULL, NULL, lb [i], ub [i], 0);

   // add auxiliary variables (unbounded for now)
   for (register int i=numcols; i<numcolsconv; i++)
     si.addCol (0, NULL, NULL, -COIN_DBL_MAX, COIN_DBL_MAX, 0);

   // get initial relaxation
   OsiCuts cs;
   couenneCg.generateCuts (si, cs);

   // store all (original + auxiliary) bounds in the relaxation
   CouNumber * colLower = new CouNumber [numcolsconv];
   CouNumber * colUpper = new CouNumber [numcolsconv];

   CouNumber *ll = couenneCg.Problem () -> Lb ();
   CouNumber *uu = couenneCg.Problem () -> Ub ();

   // overwrite original bounds, could be improved within generateCuts
   for (register int i = numcolsconv; i--;) {
     colLower [i] = (ll [i] > - COUENNE_INFINITY) ? ll [i] : -COIN_DBL_MAX;
     colUpper [i] = (uu [i] <   COUENNE_INFINITY) ? uu [i] :  COIN_DBL_MAX;
   }

   int numrowsconv = cs.sizeRowCuts ();

   // create matrix and other stuff
   CoinBigIndex * start = new CoinBigIndex [numrowsconv + 1];

   int    * length   = new int    [numrowsconv];
   double * rowLower = new double [numrowsconv];
   double * rowUpper = new double [numrowsconv];

   start[0] = 0;
   int nnz = 0;
   /* fill the four arrays. */
   for(int i = 0 ; i < numrowsconv ; i++)
   {
     OsiRowCut * cut = cs.rowCutPtr (i);

     const CoinPackedVector &v = cut->row();
     start[i+1] = start[i] + v.getNumElements();
     nnz += v.getNumElements();
     length[i] = v.getNumElements();

     rowLower[i] = cut->lb();
     rowUpper[i] = cut->ub();
   }
   assert(nnz == start[numrowsconv]);
   /* Now fill the elements arrays. */
   int * ind = new int[start[numrowsconv]];
   double * elem = new double[start[numrowsconv]];
   for(int i = 0 ; i < numrowsconv ; i++) {

     OsiRowCut * cut = cs.rowCutPtr (i);

     const CoinPackedVector &v = cut->row();

     if(v.getNumElements() != length[i])
       std::cout<<"Empty row"<<std::endl;
     //     cut->print();
     CoinCopyN(v.getIndices(), length[i], ind + start[i]);
     CoinCopyN(v.getElements(), length[i], elem + start[i]);
   }

   // Ok everything done now create interface
   CoinPackedMatrix A;
   A.assignMatrix(false, numcolsconv, numrowsconv,
                  start[numrowsconv], elem, ind,
                  start, length);
   if(A.getNumCols() != numcolsconv || A.getNumRows() != numrowsconv){
     std::cout<<"Error in row number"<<std::endl;
   }
   assert(A.getNumElements() == nnz);
   // Objective function
   double * obj = new double[numcolsconv];
   CoinFillN(obj,numcolsconv,0.);

   // some instances have no (or null) objective function, check it here
   if (couenneCg. Problem () -> nObjs () > 0)
     couenneCg.Problem() -> fillObjCoeff (obj);

   // Finally, load interface si with the initial LP relaxation
   si.loadProblem (A, colLower, colUpper, obj, rowLower, rowUpper);

   delete [] rowLower; 
   delete [] rowUpper;
   delete [] colLower;
   delete [] colUpper;
   delete [] obj;
   //   delete [] x0;

   for(int i = 0 ; i < numcols ; i++)
   {
     if(isInteger(i)){
       si.setInteger(i);
     }
   }
 
   //si.writeMpsNative("toto",NULL,NULL,1);
   si.writeLp ("toto");
   app_ -> enableWarmStart();

   //   if (problem () -> x_sol ()) {
   setColSolution (problem () -> x_sol     ());
   setRowPrice    (problem () -> duals_sol ());
     //   }
}


/** To set some application specific defaults. */
void CouenneInterface::setAppDefaultOptions(Ipopt::SmartPtr<Ipopt::OptionsList> Options){
  Options->SetStringValue("bonmin.algorithm", "B-Couenne", true, true);
  Options->SetIntegerValue("bonmin.filmint_ecp_cuts", 1, true, true);
}
} /** End Bonmin namespace. */

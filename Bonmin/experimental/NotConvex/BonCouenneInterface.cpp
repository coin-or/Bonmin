// (C) Copyright International Business Machines Corporation (IBM and Carnegie Mellon University 2006 
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 02/19/2006



#include "BonCouenneInterface.hpp"
#include "CoinHelperFunctions.hpp"


namespace Bonmin{
/** Default constructor. */
CouenneInterface::CouenneInterface():
  AmplInterface(),
  couenneCg_(NULL)
{}

/** Constructor with inputed ampl command line.*/
CouenneInterface::CouenneInterface(char **& amplArgs, SmartPtr<TNLPSolver> app):
  AmplInterface(amplArgs, app),
  couenneCg_(NULL)
  {
    const ASL_pfgh* asl = amplModel()->AmplSolverObject();
    ASL_pfgh * nc_asl = const_cast< ASL_pfgh *>(asl);
    couenneCg_ = new CouenneCutGenerator 
                       (nc_asl, false, CURRENT_ONLY,1);
  }

/** Copy constructor. */
CouenneInterface::CouenneInterface(const CouenneInterface &other):
  AmplInterface(other),
  couenneCg_(NULL)
  {
    const ASL_pfgh* asl = amplModel()->AmplSolverObject();
    ASL_pfgh * nc_asl = const_cast< ASL_pfgh *>(asl);
    couenneCg_ = new CouenneCutGenerator 
                       (nc_asl, false, CURRENT_ONLY,1);
  }

/** virutal copy constructor. */
CouenneInterface * CouenneInterface::clone(bool CopyData){
  return new CouenneInterface(*this);
}

/** Destructor. */
CouenneInterface::~CouenneInterface(){
  if(couenneCg_) delete couenneCg_;
}

/** \name Overloaded methods to build outer approximations */
  //@{
  /** \brief Extract a linear relaxation of the MINLP.
   * Solve the continuous relaxation and takes first-order outer-approximation constraints at the optimum.
   * The put everything in an OsiSolverInterface.
   */
void 
CouenneInterface::extractLinearRelaxation(OsiSolverInterface &si, bool getObj, bool solveNlp)
{
  if(solveNlp)
    initialSolve();

   // Check that couenneCg_ has been built. 
   if(couenneCg_ == NULL){
     throw CoinError("No couenne generator has been built, probably ampl .nl file was not properly read","extractLinearRelaxation", "Bonmin::CouenneInterface");
   }
   int numcols = getNumCols();
   int numcolsconv = couenneCg_->getnvars();

   CouNumber * x0 = new CouNumber[numcolsconv];
   CouNumber * colLower = new CouNumber[numcolsconv];
   CouNumber * colUpper = new CouNumber[numcolsconv];
   assert(numcolsconv >= numcols);
   
   CoinCopyN(couenneCg_->X(), numcolsconv, x0);
   CoinCopyN(couenneCg_->Lb(), numcolsconv, colLower);
   CoinCopyN(couenneCg_->Ub(), numcolsconv, colUpper);

   // Feed in the initial relaxation point
   CoinCopyN(getColSolution(), numcols, x0);
   couenneCg_->updateConv(x0, colLower, colUpper);

   int numrowsconv = couenneCg_->getncuts();

   /* Now, create matrix and other stuff. */
   CoinBigIndex * start = new CoinBigIndex[numrowsconv + 1];
   int * length = new int[numrowsconv];
   double * rowLower = new double[numrowsconv];
   double * rowUpper = new double[numrowsconv];

   start[0] = 0;
   /* fill the four arrays. */
   for(int i = 0 ; i < numrowsconv ; i++)
   {
     OsiRowCut * cut = couenneCg_->getCut(i);

     const CoinPackedVector &v = cut->row();
     start[i+1] = start[i] + v.getNumElements();
     length[i] = v.getNumElements();
     rowLower[i] = cut->lb();
     rowUpper[i] = cut->ub();
   }
   
   /* Now fill the elements arrays. */
   int * ind = new int[start[numrowsconv]];
   double * elem = new double[start[numrowsconv]];
   for(int i = 0 ; i < numrowsconv ; i++)
   {
     OsiRowCut * cut = couenneCg_->getCut(i);
     const CoinPackedVector &v = cut->row();
     CoinCopyN(v.getIndices(), length[i], ind + start[i]);
     CoinCopyN(v.getElements(), length[i], elem + start[i]);
   }

   // Ok everything done now create interface
   CoinPackedMatrix A;
   A.assignMatrix(true, numrowsconv, numcolsconv,
                  start[numrowsconv], elem, ind,
                  start, length);

   // Objective function
   double * obj = new double[numcolsconv];
   CoinFillN(obj,numcolsconv,0.);
   obj[couenneCg_->Problem()->Obj(0)->Body()->Index()] = 1;
  
   si.loadProblem(A, colLower, colUpper, obj, rowLower, rowUpper);
  

   delete [] rowLower; 
   delete [] rowUpper;
   delete [] colLower;
   delete [] colUpper;
   delete [] obj;

   for(int i = 0 ; i < numcols ; i++)
   {
     if(isInteger(i)){
	 std::cout<<"Adding integer variable :"<<i<<std::endl;
       si.setInteger(i);
     }
   }
  app_->enableWarmStart();
  setColSolution(problem()->x_sol());
  setRowPrice(problem()->duals_sol());
}

/** Get the outer approximation constraints at the currently stored optimal point.
   (Only get outer-approximations of nonlinear constraints of the problem.)*/
void 
CouenneInterface::getOuterApproximation(OsiCuts &cs, bool getObj){
   // Check that couenneCg_ has been built. 
   if(couenneCg_ == NULL){
     throw CoinError("No couenne generator has been built,"
                     " probably ampl .nl file was not properly"
                     " read",
                     "extractLinearRelaxation", 
                     "Bonmin::CouenneInterface");
   }
   int numcols = getNumCols();
   int numcolsconv = couenneCg_->getnvars();
  
   CouNumber * x = new CouNumber[numcolsconv];
   CouNumber * colLower = new CouNumber[numcolsconv];
   CouNumber * colUpper = new CouNumber[numcolsconv];
   assert(numcolsconv >= numcols);
   
   CoinCopyN(couenneCg_->X(), numcolsconv, x);
   CoinCopyN(couenneCg_->Lb(), numcolsconv, colLower);
   CoinCopyN(couenneCg_->Ub(), numcolsconv, colUpper);

   // Feed in the current state
   CoinCopyN(getColSolution(), numcols, x);
   CoinCopyN(getColLower(), numcols, colLower);
   CoinCopyN(getColUpper(), numcols, colUpper);

   couenneCg_->updateConv(x, colLower, colUpper);


   int ncuts = couenneCg_->getncuts();
   for(int i = 0 ; i < ncuts ; i++)
   {
     cs.insert(*couenneCg_->getCut(i)); 
   }
   
   delete [] colLower;
   delete [] colUpper;
   delete [] x;
}

} /** End Bonmin namespace. */

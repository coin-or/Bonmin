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
#include "CouenneProblem.hpp"

namespace Bonmin {

/** Default constructor. */
CouenneInterface::CouenneInterface():
  AmplInterface(),
  have_nlp_solution_ (false)
{}

/** Copy constructor. */
CouenneInterface::CouenneInterface(const CouenneInterface &other):
  AmplInterface(other),
  have_nlp_solution_ (false)
{}

/** virtual copy constructor. */
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

  CouenneProblem *p = couenneCg.Problem ();
  bool is_feasible = true;

  if (solveNlp) {

    int nvars = p -> nVars();

    if (p -> doFBBT ()) {

      // include the rhs of auxiliary-based constraints into the FBBT
      // (should be useful with Vielma's problems, for example)

      for (int i=0; i < p -> nCons (); i++) {

	// for each constraint
	CouenneConstraint *con = p -> Con (i);

	// (which has an aux as its body)
	int index = con -> Body () -> Index ();

	if ((index >= 0) && (con -> Body () -> Type () == AUX)) {

	  // if there exists violation, add constraint
	  CouNumber 
	    l = con -> Lb () -> Value (),	
	    u = con -> Ub () -> Value ();

	  // tighten bounds in Couenne's problem representation
	  p -> Lb (index) = CoinMax (l, p -> Lb (index));
	  p -> Ub (index) = CoinMin (u, p -> Ub (index));
	}
      }

      t_chg_bounds *chg_bds = new t_chg_bounds [nvars];

      for (int i=0; i<nvars; i++) {
	chg_bds [i].setLower(t_chg_bounds::CHANGED);
	chg_bds [i].setUpper(t_chg_bounds::CHANGED);
      }

      if (!(p -> boundTightening (chg_bds, NULL))) {
	is_feasible = false;
	printf ("warning, tightened NLP is infeasible\n");
      }

      delete [] chg_bds;

      const double 
	*nlb = getColLower (),
	*nub = getColUpper ();

      for (int i=0; i < p -> nOrigVars (); i++) 
	if (p -> Var (i) -> Multiplicity () > 0) {
	  if (nlb [i] < p -> Lb (i) - COUENNE_EPS) setColLower (i, p -> Lb (i));
	  if (nub [i] > p -> Ub (i) + COUENNE_EPS) setColUpper (i, p -> Ub (i));
	} else { 
	  // if not enabled, fix them in the NLP solver
	  setColLower (i, 0.);
	  setColUpper (i, 0.);
	}
    }

    if (is_feasible) 
      initialSolve ();
    else {
      OsiAuxInfo * auxInfo = si.getAuxiliaryInfo ();
      BabInfo * babInfo = dynamic_cast <BabInfo *> (auxInfo);

      if (babInfo) 
	babInfo -> setInfeasibleNode ();
    }

    if (is_feasible && isProvenOptimal ()) {

      CouNumber obj             = getObjValue    ();
      const CouNumber *solution = getColSolution ();

      if (getNumIntegers () > 0) {

	int
	  norig = p -> nOrigVars (),
	  nvars = p -> nVars ();

	bool fractional = false;

	// if problem is integer, check if any integral variable is
	// fractional. If so, round them and re-optimize

	for (int i=0; i<norig; i++)
	  if (p  -> Var (i) -> isInteger () &&
	      (p -> Var (i) -> Multiplicity () > 0) &&
	      (!::isInteger (solution [i]))) {
	    fractional = true;
	    break;
	  }

	if (fractional) { // try again if solution found by Ipopt is fractional

	  double 
	    *lbSave = new double [norig],
	    *ubSave = new double [norig],

	    *lbCur  = new double [nvars],
	    *ubCur  = new double [nvars],

	    *Y      = new double [nvars];

	  CoinCopyN (getColLower (), norig, lbSave);
	  CoinCopyN (getColUpper (), norig, ubSave);

	  CoinFillN (Y,     nvars, 0.);
	  CoinFillN (lbCur, nvars, -COUENNE_INFINITY);
	  CoinFillN (ubCur, nvars,  COUENNE_INFINITY);

	  CoinCopyN (getColLower (), norig, lbCur);
	  CoinCopyN (getColUpper (), norig, ubCur);

	  if (p -> getIntegerCandidate (solution, Y, lbCur, ubCur) >= 0) {

	    for (int i = getNumCols (); i--;) {

	      if (lbCur [i] > ubCur [i]) {
		double swap = lbCur [i];
		lbCur [i] = ubCur [i];
		ubCur [i] = swap;
	      }

	      if      (Y [i] < lbCur [i]) Y [i] = lbCur [i];
	      else if (Y [i] > ubCur [i]) Y [i] = ubCur [i];
	    }

	    for (int i=0; i<norig; i++)
	      if ((p -> Var (i) -> Multiplicity () > 0) &&
		  p  -> Var (i) -> isInteger ()) {
		setColLower (i, lbCur [i]);
		setColUpper (i, ubCur [i]);
	      }

	    setColSolution (Y); // use initial solution given 

	    resolve (); // solve with integer variables fixed

	    obj      = getObjValue ();
	    solution = getColSolution ();

	    // restore previous bounds on integer variables
	    for (int i=0; i<norig; i++)
	      if ((p -> Var (i) -> Multiplicity () > 0) &&
		  p  -> Var (i) -> isInteger ()) {
		setColLower (i, lbSave [i]);
		setColUpper (i, ubSave [i]);
	      }
	  }

	  delete [] Y;
	  delete [] lbSave;
	  delete [] ubSave;
	  delete [] lbCur;
	  delete [] ubCur;
	} 
      }

      // re-check optimality in case resolve () was called
      if (isProvenOptimal () && 
	  p -> checkNLP (solution, obj, true) && // true for recomputing obj
	  (obj < p -> getCutOff ())) {

	// tell caller there is an initial solution to be fed to the initHeuristic
	have_nlp_solution_ = true;

	// set cutoff to take advantage of bound tightening
	p -> setCutOff (obj);

	OsiAuxInfo * auxInfo = si.getAuxiliaryInfo ();
	BabInfo * babInfo = dynamic_cast <BabInfo *> (auxInfo);

	if (babInfo) {
	  babInfo -> setNlpSolution (solution, getNumCols (), obj);
	  babInfo -> setHasNlpSolution (true);
	}
      }
    }
  }

  if (!is_feasible) // nothing else to do, problem infeasible
    return;

  int 
    numcols     = getNumCols (), // # original               variables
    numcolsconv = p -> nVars (); // # original + # auxiliary variables

  const double
    *lb = getColLower (),
    *ub = getColUpper ();

   // add original and auxiliary variables to the new problem
   for (int i=0; i<numcols; i++) 
     if (p -> Var (i) -> Multiplicity () > 0) si.addCol (0, NULL,NULL, lb [i],       ub [i],      0);
     else                                     si.addCol (0, NULL,NULL, 0.,           0.,          0);
   for (int i=numcols; i<numcolsconv; i++)    si.addCol (0, NULL,NULL, -COIN_DBL_MAX,COIN_DBL_MAX,0);

   // get initial relaxation
   OsiCuts cs;
   couenneCg.generateCuts (si, cs);

   // store all (original + auxiliary) bounds in the relaxation
   CouNumber * colLower = new CouNumber [numcolsconv];
   CouNumber * colUpper = new CouNumber [numcolsconv];

   CouNumber *ll = p -> Lb ();
   CouNumber *uu = p -> Ub ();

   // overwrite original bounds, could have improved within generateCuts
   for (int i = numcolsconv; i--;) 
     if (p -> Var (i) -> Multiplicity () > 0) {
       colLower [i] = (ll [i] > - COUENNE_INFINITY) ? ll [i] : -COIN_DBL_MAX;
       colUpper [i] = (uu [i] <   COUENNE_INFINITY) ? uu [i] :  COIN_DBL_MAX;
     } else {
       colLower [i] = 0.;
       colUpper [i] = 0.;
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

   assert (nnz == start [numrowsconv]);
   /* Now fill the elements arrays. */
   int * ind = new int[start[numrowsconv]];
   double * elem = new double[start[numrowsconv]];
   for(int i = 0 ; i < numrowsconv ; i++) {

     OsiRowCut * cut = cs.rowCutPtr (i);

     const CoinPackedVector &v = cut->row();

     if(v.getNumElements() != length[i])
       std::cout<<"Empty row"<<std::endl;
     //     cut->print();
     CoinCopyN (v.getIndices(),  length[i], ind  + start[i]);
     CoinCopyN (v.getElements(), length[i], elem + start[i]);
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
   if (p -> nObjs () > 0)
     p -> fillObjCoeff (obj);

   // Finally, load interface si with the initial LP relaxation
   si.loadProblem (A, colLower, colUpper, obj, rowLower, rowUpper);

   delete [] rowLower; 
   delete [] rowUpper;
   delete [] colLower;
   delete [] colUpper;
   delete [] obj;

   for (int i=0; i<numcolsconv; i++)
     if ((p -> Var (i) -> Multiplicity () > 0) &&
	 (p -> Var (i) -> isDefinedInteger ()))
       si.setInteger (i);
 
   //si.writeMpsNative("toto",NULL,NULL,1);
   //si.writeLp ("toto");
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

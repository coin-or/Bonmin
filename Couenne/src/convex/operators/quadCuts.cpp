/*
 * Name:    quadCuts.cpp
 * Authors: Pierre Bonami
 *          Stefan Vigerske
 *          Pietro Belotti
 * Purpose: based on upper and lower convexification, add cuts to convexify
 *
 * (C) International Business Machines 2007. This file is licensed under the Common Public License (CPL)
 */

#include <exprQuad.hpp>

#include <CouenneProblem.hpp>
#include <CouenneCutGenerator.hpp>

#include "CoinHelperFunctions.hpp"
//#define DEBUG
void exprQuad::quadCuts(exprAux *w, OsiCuts &cs, const CouenneCutGenerator *cg){

  assert(dIndex_ != NULL);

#ifdef DEBUG
  std::cout<<"Expression has "<<nlterms_<<" linear terms and "
           <<nqterms_<<" quadratic terms."<<std::endl;
#endif

  //First get on which side constraint is violated to get the good lambda
  double * lambda = NULL;
  double exprVal = (*this) ();
  double varVal =  expression::Variable (w -> Index());
  if(varVal < exprVal)//Use under-estimator
  {
     lambda = dCoeffLo_;
  }
  else //use overestimator
  {
     lambda = dCoeffUp_;
  }
  const CouenneProblem& problem = *(cg->Problem());
  const int & numcols = problem.nOrig();
  const double * colsol = problem.X();
  const double * lower = problem.Lb();
  const double * upper = problem.Ub();

  std::cout<<"Point to cut"<<std::endl;
  for(int i = 0 ; i < numcols ; i++){
    std::cout<<colsol[i]<<", ";}
  std::cout<<std::endl;

  //Initialize by copying a into a dense vector and computing Q x^*
  double * vec = new double[numcols];
  CoinFillN(vec, numcols, 0.);

  //Start by computing Q x^*.

  for(int k = 0 ; k < nqterms_ ; k++) {
     vec[qindexI_[k]] += qcoeff_[k] * colsol[qindexJ_[k]];
     vec[qindexJ_[k]] += qcoeff_[k] * colsol[qindexI_[k]];
  }

  //multiply it by x^*^T again and store the result for the lower bound
  double a0 = - c0_;
  for(int i = 0 ; i < numcols ; i++){
    a0 += vec[i] * colsol[i];
    vec[i] *= 2;
  }

  // Add a to it.
  for(int i = 0 ; i < nlterms_ ; i++){
     vec[index_[i]] += coeff_[i];
  }


  // Don't forget the auxiliaries!
  /*for(int i = 0 ; i < nargs_ ; i++){
     vec[arglist_ [i] -> Index ()] += 1;
     }*/

  // And myself
  vec [w -> Index ()] -= 1;

  if (lambda != NULL) {
    // Now the part which depends on lambda
    for(int k = 0 ; k < nDiag_ ; k++){
      a0 += lambda[k] * lower[dIndex_[k]] * upper[dIndex_[k]];
      a0 += lambda[k] * colsol[dIndex_[k]] * colsol[dIndex_[k]];
      vec[dIndex_[k]] += lambda[k] * (lower[dIndex_[k]] + upper[dIndex_[k]]);
      vec[dIndex_[k]] -= lambda[k] * (colsol[dIndex_[k]]) * 2;
    }
  }



  // Count the number of non-zeroes
  int nnz = 0;
  for(int i = 0 ; i < numcols ; i++){
    if(fabs(vec[i]) > COUENNE_EPS){
       nnz++;
    }
  }
#ifdef DEBUG
  std::cout<<"My cut should have "<<nnz<<" non zeroes."<<std::endl;
#endif
  // Pack the vector into a CoinPackedVector and generate the cut.
  CoinPackedVector a(false);
  a.reserve(nnz);
  for(int i = 0 ; i < numcols ; i++){
    if(fabs(vec[i]) > COUENNE_EPS){
       a.insert(i,vec[i]);
    }
  }
  OsiRowCut cut;
  cut.setRow(a);
  if( lambda == dCoeffLo_){
     cut.setUb(a0);
     cut.setLb(-COUENNE_INFINITY);
  }
  else {
    cut.setLb(a0);
    cut.setUb(COUENNE_INFINITY);
  }
  delete [] vec;
  cs.insert(cut);
}

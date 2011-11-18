// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 10/16/2007
#include "BonTMINLP2OsiLP.hpp"
#include "BonTypes.hpp"
#include "OsiSolverInterface.hpp"
#include "BonTMINLP2TNLP.hpp"
#include "CoinPackedMatrix.hpp"

#include <vector>
#include <sstream>
#include <climits>

using namespace Ipopt;

namespace Bonmin {

int TMINLP2OsiLP::nTimesCalled = 0;

void
TMINLP2OsiLP::extract(OsiSolverInterface *si, 
                     const double * x, bool getObj)
{
  assert(IsValid(model_));
  int n;
  int m;
  int nnz_jac_g;
  int nnz_h_lag;
  TNLP::IndexStyleEnum index_style;
  //Get problem information
  model_->get_nlp_info( n, m, nnz_jac_g, nnz_h_lag, index_style);

  //get Jacobian
  model_->eval_jac_g(n, x, 1, m, nnz_jac_g, NULL, NULL, value_());


  vector<double> g(m);
  model_->eval_g(n, x, 1, m, g());

  vector<double> rowLow(m);
  vector<double> rowUp(m);



  const double * rowLower = model_->g_l();
  const double * rowUpper = model_->g_u();
  const double * colLower = model_->x_l();
  const double * colUpper = model_->x_u();

  double infty = si->getInfinity();
  double nlp_infty = infty_;

  for(int i = 0 ; i < m ; i++) {
    if(const_types_[i] == Ipopt::TNLP::NON_LINEAR) {
      if(rowLower[i] > - nlp_infty){
        rowLow[i] = (rowLower[i] - g[i]) - 1e-07;
      }
      else
        rowLow[i] = - infty;
      if(rowUpper[i] < nlp_infty)
        rowUp[i] =  (rowUpper[i] - g[i]) + 1e-07;
      else
        rowUp[i] = infty;
    }
    else {
      if(rowLower[i] > -nlp_infty){
         rowLow[i] = (rowLower[i]);
      }
      else
        rowLow[i] = - infty;
      if(rowUpper[i] < nlp_infty){
         rowUp[i] =  (rowUpper[i]);
      }
      else
        rowUp[i] = infty;
    }
  }

  
  
  //Then convert everything to a CoinPackedMatrix
  //Go through values, clean coefficients and fix bounds
  for(int i = 0 ; i < nnz_jac_g ; i++) {
    if(const_types_[iRow_[i]] != TNLP::LINEAR){//For linear just copy is fine.
       if(//For other clean tinys
       cleanNnz(value_[i],colLower[jCol_[i]], colUpper[jCol_[i]],
                rowLower[iRow_[i]], rowUpper[iRow_[i]],
                x[jCol_[i]],
                rowLow[iRow_[i]],
                rowUp[iRow_[i]], tiny_, very_tiny_)) {      
          if(rowLow[iRow_[i]] > -infty)
            rowLow[iRow_[i]] += value_[i] * x[jCol_[i]];
          if(rowUp[iRow_[i]] < infty)
            rowUp[iRow_[i]] += value_[i] *x[jCol_[i]];
       }
    }
  }
  CoinPackedMatrix mat(true, iRow_(), jCol_(), value_(), nnz_jac_g);
  mat.setDimensions(m,n); // In case matrix was empty, this should be enough
  
#if 0
  vector<double> act(m);
  mat.times(x, act());
  for(int j = 0 ; j < m ; j++){
    if(j==10  && fabs(x[0] - 4.73032) < 1e-4)
    assert(act[j] + 1e-5 > rowLow[j] && act[j] - 1e-5 < rowUp[j] + g[j]);
  }
#endif

  vector<double> obj(n);
  for(int i = 0 ; i < n; i++)
    obj[i] = 0.;
  
  
  si->loadProblem(mat, colLower, colUpper, obj(), rowLow(), rowUp());
  for(int i = 0 ; i < n ; i++) {
    if(model_->var_types()[i] == TMINLP::BINARY || model_->var_types()[i] == TMINLP::INTEGER )
      si->setInteger(i);
  }
  if(getObj) {
     bool addObjVar = false;
     if(model_->hasLinearObjective()){
       double zero;
       vector<double> x0(n,0.);
       model_->eval_f(n, x0(), 1, zero);
       si->setDblParam(OsiObjOffset, -zero);
       //Copy the linear objective and don't create a dummy variable.
       model_->eval_grad_f(n, x, 1,obj());
       si->setObjective(obj());
    }
    else {
      addObjVar = true;
   }
   
   if(addObjVar){
      addObjectiveFunction(*si, x);
    }
  }

}
 
void TMINLP2OsiLP::get_oas(OsiCuts &cs, const double *x, bool getObj, bool global) {

  int n,m, nnz_jac_g, nnz_h_lag;
  TNLP::IndexStyleEnum index_style;
  model_->get_nlp_info( n, m, nnz_jac_g, nnz_h_lag, index_style);


  vector<double> g(m);

  model_->eval_jac_g(n, x, 1, m, nnz_jac_g, NULL, NULL, value_());
  model_->eval_g(n,x,1,m,g());

  //As jacobian is stored by cols fill OsiCuts with cuts
  vector<double> lb(m + 1);
  vector<double> ub(m + 1);

  vector<int> row2cutIdx(m,-1);//store correspondance between index of row and index of cut (some cuts are not generated for rows because linear, or not binding). -1 if constraint does not generate a cut, otherwise index in cuts.
  int numCuts = 0;

  const double * rowLower = model_->g_l();
  const double * rowUpper = model_->g_u();
  const double * colLower = model_->x_l();
  const double * colUpper = model_->x_u();
  double infty = infty_;
  double nlp_infty = infty_;
  
  for(int rowIdx = 0; rowIdx < m ; rowIdx++) {
    if(const_types_[rowIdx] == TNLP::NON_LINEAR) {
      row2cutIdx[rowIdx] = numCuts;
      if(rowLower[rowIdx] > - nlp_infty)
        lb[numCuts] = rowLower[rowIdx] - g[rowIdx];
      else
        lb[numCuts] = - infty;
      if(rowUpper[rowIdx] < nlp_infty)
        ub[numCuts] = rowUpper[rowIdx] - g[rowIdx];
      else
        ub[numCuts] = infty;
      numCuts++;
    }
  }

  lb.resize(numCuts);
  ub.resize(numCuts);
  vector<CoinPackedVector>  cuts(numCuts);


  for(int i = 0 ; i < nnz_jac_g ; i++) {
    const int &rowIdx = iRow_[i];
    const int & cutIdx = row2cutIdx[ rowIdx ];
    if(cutIdx != -1) {
      const int &colIdx = jCol_[i];
      //"clean" coefficient
      if(cleanNnz(value_[i],colLower[colIdx], colUpper[colIdx],
		  rowLower[rowIdx], rowUpper[rowIdx],
		  x[colIdx],
		  lb[cutIdx],
		  ub[cutIdx], tiny_, very_tiny_)) {
        cuts[cutIdx].insert(colIdx,value_[i]);
        if(lb[cutIdx] > - infty)
          lb[cutIdx] += value_[i] * x[colIdx];
        if(ub[cutIdx] < infty)
	  ub[cutIdx] += value_[i] * x[colIdx];
      }
    }
  }

  for(int cutIdx = 0; cutIdx < numCuts; cutIdx++) {
    OsiRowCut newCut;
    if(global)
      newCut.setGloballyValidAsInteger(1);
    newCut.setLb(lb[cutIdx]);
    newCut.setUb(ub[cutIdx]);
    newCut.setRow(cuts[cutIdx]);
    cs.insert(newCut);
  }

  if(getObj && !model_->hasLinearObjective()) { // Get the objective cuts
    vector<double> obj(n);
    model_->eval_grad_f(n, x, 1,obj());
    double f;
    model_->eval_f(n, x, 1, f);

    CoinPackedVector v;
    v.reserve(n);
    lb.back() = -f;
    ub.back() = -f;
    //double minCoeff = 1e50;
    for(int i = 0; i<n ; i++) {
      if(cleanNnz(obj[i],colLower[i], colUpper[i],
          -infty_, 0,
          x[i],
          lb.back(),
          ub.back(),tiny_, very_tiny_)) {
        //	      minCoeff = min(fabs(obj[i]), minCoeff);
        v.insert(i,obj[i]);
        lb.back() += obj[i] * x[i];
        ub.back() += obj[i] * x[i];
      }
    }
    v.insert(n,-1);
    //Compute cut violation
    OsiRowCut newCut;
    newCut.setGloballyValidAsInteger(1);
    newCut.setRow(v);
    newCut.setLb(-DBL_MAX);
    newCut.setUb(ub.back());
    cs.insert(newCut);
    }
    }

void 
TMINLP2OsiLP::addObjectiveFunction(OsiSolverInterface &si, 
                                const double *x){
  assert(IsValid(model_));
  const double * colLower = model_->x_l();
  const double * colUpper = model_->x_u();
  int numcols = model_->num_variables();
  assert(numcols == si.getNumCols() );
  vector<double> obj(numcols);
  model_->eval_grad_f(numcols, x, 1,obj());
  //add variable alpha
  //(alpha should be empty in the matrix with a coefficient of -1 and unbounded)
  CoinPackedVector a;
  si.addCol(a,-si.getInfinity(), si.getInfinity(), 1.);
  
  // Now get the objective cuts
  // get the gradient, pack it and add the cut
  double ub;
  model_->eval_f(numcols, x, 1, ub);
  ub*=-1;
  double lb = -infty_;
  CoinPackedVector objCut;
  CoinPackedVector * v = &objCut;
  v->reserve(numcols+1);
  for(int i = 0; i<numcols ; i++) {
   if(si.getNumRows())
   {
    if(cleanNnz(obj[i],colLower[i], colUpper[i],
        -infty_, 0,
        x[i],
        lb,
        ub, tiny_, very_tiny_)) {
      v->insert(i,obj[i]);
      lb += obj[i] * x[i];
      ub += obj[i] * x[i];
    }
   }
   else //Unconstrained problem can not put clean coefficient
   {
       if(cleanNnz(obj[i],colLower[i], colUpper[i],
        -infty_, 0,
        x[i],
        lb,
        ub, 1e-03, 1e-08)) {
      v->insert(i,obj[i]);
      lb += obj[i] * x[i];
      ub += obj[i] * x[i];
       }
   }
  }
  v->insert(numcols,-1);
  si.addRow(objCut, lb, ub);
  }


   void 
   TMINLP2OsiLP::initialize_jac_storage(){
     assert(IsValid(model_));
     int n;
     int m;
     int nnz_jac_g;
     int nnz_h_lag;
     TNLP::IndexStyleEnum index_style;
     //Get problem information
     model_->get_nlp_info( n, m, nnz_jac_g, nnz_h_lag, index_style);
     jCol_.resize(nnz_jac_g);
     iRow_.resize(nnz_jac_g);
     value_.resize(nnz_jac_g);

     model_->eval_jac_g(n, NULL, 0, m, nnz_jac_g, iRow_(), jCol_(), NULL);
     if(index_style == TNLP::FORTRAN_STYLE){
       for(size_t i = 0 ; i < iRow_.size() ; i++){
         iRow_[i]--; jCol_[i]--;
       }
     }

     const_types_.resize(m);
     model_->get_constraints_linearity(m, const_types_());

   }

}


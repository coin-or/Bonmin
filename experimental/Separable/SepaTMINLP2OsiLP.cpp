// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 10/16/2007
#include "SepaTMINLP2OsiLP.hpp"
#include "BonTypes.hpp"
#include "OsiSolverInterface.hpp"
#include "BonTMINLP2TNLP.hpp"
#include "CoinPackedMatrix.hpp"

#include <vector>
#include <sstream>
#include <climits>

using namespace Ipopt;

namespace Sepa {

void
SepaTMINLP2OsiLP::extract(OsiSolverInterface *si, 
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


  Bonmin::vector<double> g(m);
  model_->eval_g(n, x, 1, m, g());

  Bonmin::vector<double> rowLow(m);
  Bonmin::vector<double> rowUp(m);



  const double * rowLower = model_->g_l();
  const double * rowUpper = model_->g_u();
  const double * colLower = model_->x_l();
  const double * colUpper = model_->x_u();

  double nlp_infty = si->getInfinity();
  double infty = DBL_MAX;

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

  Bonmin::vector<double> obj(n);
  for(int i = 0 ; i < n; i++)
    obj[i] = 0.;
  
  
  si->loadProblem(mat, colLower, colUpper, obj(), rowLow(), rowUp());
  for(int i = 0 ; i < n ; i++) {
    if(model_->var_types()[i] == Bonmin::TMINLP::BINARY || model_->var_types()[i] == Bonmin::TMINLP::INTEGER )
      si->setInteger(i);
  }
  if(getObj) {
     if(model_->hasLinearObjective()){
       double zero;
       Bonmin::vector<double> x0(n,0.);
       model_->eval_f(n, x0(), 1, zero);
       si->setDblParam(OsiObjOffset, -zero);
       //Copy the linear objective and don't create a dummy variable.
       model_->eval_grad_f(n, x, 1,obj());
       si->setObjective(obj());
    }
    else {
      throw -1;
   }
   
  }
  
  OsiCuts cs;
  get_oas(cs, x, 0, 1);
  si->applyCuts(cs);

  
}
 
void 
SepaTMINLP2OsiLP::get_oas(OsiCuts &cs, const double *x, bool getObj, bool global) const {

  int n,m, nnz_jac_g, nnz_h_lag;
  TNLP::IndexStyleEnum index_style;
  model_->get_nlp_info( n, m, nnz_jac_g, nnz_h_lag, index_style);

  Bonmin::vector<double> g(m);

  model_->eval_jac_g(n, x, 1, m, nnz_jac_g, NULL, NULL, value_());
  model_->eval_g(n,x,0,m,g());

  //As jacobian is stored by cols fill OsiCuts with cuts
  Bonmin::vector<double> lb(m + 1);
  Bonmin::vector<double> ub(m + 1);

  Bonmin::vector<int> row2cutIdx(m,-1);//store correspondance between index of row and index of cut (some cuts are not generated for rows because linear, or not binding). -1 if constraint does not generate a cut, otherwise index in cuts.

  std::vector<int> cut2rowIdx;

  int numCuts = 0;

  const double * rowLower = model_->g_l();
  const double * rowUpper = model_->g_u();
  const double * colLower = model_->x_l();
  const double * colUpper = model_->x_u();
  double nlp_infty = infty_;
  double infty = DBL_MAX;
  
  for(int rowIdx = 0; rowIdx < m ; rowIdx++) {
    if(const_types_[rowIdx] == TNLP::NON_LINEAR) {
      row2cutIdx[rowIdx] = numCuts;
      cut2rowIdx.push_back(rowIdx);
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
  Bonmin::vector<CoinPackedVector>  cuts(numCuts);


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
        if(fabs(value_[i]) > 1e15) {
           printf("Not generating cut because of big coefficient %g col %i x[%i] = %g\n", value_[i], colIdx, colIdx, x[colIdx]);
           return;
        }
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
    //********* Perspective Extension ********//
    const int* ids = model_->get_const_xtra_id(); // vector of indices corresponding to the binary variable activating the corresponding constraint
    // Get the index of the corresponding indicator binary variable
    int binary_id = (ids == NULL) ? -1 : ids[cut2rowIdx[cutIdx]];// index corresponding to the binary variable activating the corresponding constraint
    if(binary_id>0) {// If this hyperplane is a linearization of a disjunctive constraint, we link its righthand (or lefthand) side to the corresponding indicator binary variable
        //printf("Using perspectives\n");
        if (lb[cutIdx] > -infty) { // ∂x ≥ lb => ∂x - lb*z ≥ 0
            cuts[cutIdx].insert(binary_id, -lb[cutIdx]);
            newCut.setLb(0);
            newCut.setUb(ub[cutIdx]);
            
        }
        if (ub[cutIdx] < infty) { // ∂x ≤ ub => ∂x - ub*z ≤ 0
            cuts[cutIdx].insert(binary_id, -ub[cutIdx]);
            newCut.setLb(lb[cutIdx]);
            newCut.setUb(0);
            
        }
    }
    else {
        newCut.setLb(lb[cutIdx]);
        newCut.setUb(ub[cutIdx]);
    }
    //********* Perspective Extension ********//
    newCut.setRow(cuts[cutIdx]);
    cs.insert(newCut);
  }
  printf("++++++++   I have generated %i cuts   +++++++++\n", numCuts);

  return;


}

void 
SepaTMINLP2OsiLP::get_oa(int rowIdx, OsiCuts &cs, const double *x, bool getObj, bool global) const {

  int n,m, nnz_jac_g, nnz_h_lag;
  TNLP::IndexStyleEnum index_style;
  model_->get_nlp_info( n, m, nnz_jac_g, nnz_h_lag, index_style);

  double gi;
  model_->eval_gi(n,x,1, rowIdx,gi);
  Bonmin::vector<int> jCol(n);
  int nnz;
  model_->eval_grad_gi(n, x, 0, rowIdx, nnz, jCol(), NULL);
  model_->eval_grad_gi(n, x, 0, rowIdx, nnz, NULL, value_());


  //As jacobian is stored by cols fill OsiCuts with cuts
  double lb;
  double ub;

  const double * rowLower = model_->g_l();
  const double * rowUpper = model_->g_u();
  const double * colLower = model_->x_l();
  const double * colUpper = model_->x_u();
  double nlp_infty = infty_;
  double infty = DBL_MAX;
  
  if (rowLower[rowIdx] > -nlp_infty)
    lb = rowLower[rowIdx] - gi;
  else
    lb = -infty;
  if (rowUpper[rowIdx] < nlp_infty)
    ub = rowUpper[rowIdx] - gi;
  else
    ub = infty;

  CoinPackedVector  cut;


  for(int i = 0 ; i < nnz ; i++) {
    if(index_style == Ipopt::TNLP::FORTRAN_STYLE) jCol[i]--;
    const int &colIdx = jCol[i];
    //"clean" coefficient
    if(cleanNnz(value_[i],colLower[colIdx], colUpper[colIdx],
      rowLower[rowIdx], rowUpper[rowIdx],
      x[colIdx],
      lb,
      ub, tiny_, very_tiny_)) {
      if(fabs(value_[i]) > 1e15) {
         printf("Not generating cut because of big coefficient %g col %i x[%i] = %g\n", value_[i], colIdx, colIdx, x[colIdx]);
         return;
      }
      cut.insert(colIdx,value_[i]);
        if(lb > - infty)
          lb += value_[i] * x[colIdx];
        if(ub < infty)
           ub += value_[i] * x[colIdx];
      }
  }

  OsiRowCut newCut;
  if(global)
    newCut.setGloballyValidAsInteger(1);
  //********* Perspective Extension ********//
  const int* ids = model_->get_const_xtra_id(); // vector of indices corresponding to the binary variable activating the corresponding constraint
  // Get the index of the corresponding indicator binary variable
  int binary_id = (ids == NULL) ? -1 : ids[rowIdx];// index corresponding to the binary variable activating the corresponding constraint
    if(binary_id>0) {// If this hyperplane is a linearization of a disjunctive constraint, we link its righthand (or lefthand) side to the corresponding indicator binary variable
        //printf("Using perspectives\n");
        if (lb > -infty) { // ∂x ≥ lb => ∂x - lb*z ≥ 0
            cut.insert(binary_id, -lb);
            newCut.setLb(0);
            newCut.setUb(ub);
            
        }
        if (ub < infty) { // ∂x ≤ ub => ∂x - ub*z ≤ 0
            cut.insert(binary_id, -ub);
            newCut.setLb(lb);
            newCut.setUb(0);
        }
    }
    else {
        newCut.setLb(lb);
        newCut.setUb(ub);
    }
    //********* Perspective Extension ********//
    newCut.setRow(cut);
    cs.insert(newCut);
  return;
}

void 
SepaTMINLP2OsiLP::get_refined_oa(OsiCuts & cs) const{
   if(num_approx_ <= 0)
     return;
   int n;
   int m;
   int nnz_jac_g;
   int nnz_h_lag;
   Ipopt::TNLP::IndexStyleEnum index_style;
   //Get problem information
   model_->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);

  //const double * rowLower = model_->g_l();
  //const double * rowUpper = model_->g_u();
  const double * colLower = model_->x_l();
  const double * colUpper = model_->x_u();
  Bonmin::vector<Ipopt::TNLP::LinearityType>  varTypes(n);

   model_->get_variables_linearity(n, varTypes());
   //const Bonmin::TMINLP::VariableType* variableType = model_->var_types();
   // Hassan OA initial description

   double * p = CoinCopyOfArray(colLower, n);
   double * pp = CoinCopyOfArray(colLower, n);
   double * up = CoinCopyOfArray(colUpper, n);

   std::vector<double> step(n);


     for (int i = 0; i < n; i++) {
       if (p[i] < -1e4){
          p[i] = pp[i] = -1e4;
       }
       if (up[i] > 1e4){
          up[i] = 1e4;
       }
     } 


   //Step step
   for (int i = 0; i < n; i++) {

      if (varTypes[i] == Ipopt::TNLP::LINEAR) {
         step[i] = 0;
         p[i] = pp[i] = up[i] = 0;
      }
      else
        step[i] = (up[i] - p[i]) / (num_approx_);

    }
    get_oas(cs, p, 0, true);// Generate Tangents at current point        
    for (int j = 1; j <= num_approx_; j++) {

      for (int i = 0; i < n; i++) {
        pp[i] += step[i];
      }
      
      get_oas(cs, pp, 0, true);// Generate Tangents at current point

   }

   get_oas(cs, up, 0, true);// Generate Tangents at current point
      
   delete [] p;
   delete [] pp;
   delete [] up;
  }

void 
SepaTMINLP2OsiLP::add_outer_description_function_values(OsiSolverInterface &si) {
   int n;
   int m;
   int nnz_jac_g;
   int nnz_h_lag;
   Ipopt::TNLP::IndexStyleEnum index_style;
   //Get problem information
   model_->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);

  //const double * rowLower = model_->g_l();
  //const double * rowUpper = model_->g_u();
  const double * colLower = model_->x_l();
  const double * colUpper = model_->x_u();
   Bonmin::vector<Ipopt::TNLP::LinearityType>  varTypes(n);

   model_->get_variables_linearity(n, varTypes());
   const Bonmin::TMINLP::VariableType* variableType = model_->var_types();
   // Hassan OA initial description
   OsiCuts cs;

   double * p = CoinCopyOfArray(colLower, n);
   double * pp = CoinCopyOfArray(colLower, n);
   double * up = CoinCopyOfArray(colUpper, n);
   double * low = CoinCopyOfArray(colLower, n);

   Bonmin::vector<double> step(n);
   Bonmin::vector<int> nbG(m,2);


     for (int i = 0; i < n; i++) {
       if (low[i] < -1e4){
          low[i] = -1e4;
       }
       if (up[i] > 1e4){
          up[i] = 1e4;
       }
     } 


   for (int i = 0; i < n; i++) {

      if (varTypes[i] == Ipopt::TNLP::LINEAR) {
         step[i] = 0;
         low[i] = p[i] = pp[i] = up[i] = 0;
      }
      else
        step[i] = (up[i] - low[i]) / (num_approx_*20);

    }
    get_oas(cs, low, 0, true);// Generate Tangents at start point        
    get_oas(cs, up, 0, true);// Generate Tangents at end point        
    Bonmin::vector<double> g_p(m);
    Bonmin::vector<double> g_pp(m);
    model_->eval_g(n, low, 1, m, g_p()); // Evaluate function g at lowerbounds
    model_->eval_g(n, up, 1, m, g_pp()); // Evaluate function g at upperbounds

    for (int i = 0; (i < m); i++) { // Generate Outer-Approximation initial cuts for all nonlinear constraints
      if(const_types_[i] != Ipopt::TNLP::NON_LINEAR) continue;
      double thresh = std::abs(g_pp[i] - g_p[i])/(num_approx_-1);
      
      printf("Constraint %i threshold is %g lower val %g upper %g\n",i, thresh,
             g_p[i], g_pp[i]);
      std::copy(low, low + n, p);
      std::copy(low, low + n, pp);
      double g_p_i, g_pp_i;
      g_p_i = g_p[i];
      int n_steps = 0;
      while (nbG[i] < num_approx_) { // Iterate untill increase is sufficient 
        n_steps++;
        // Curvature sampling
        for (int j = 0; j < n; j++) {
          pp[j] += step[j];
        }
        model_->eval_gi(n, pp, 1, i, g_pp_i);
    //    printf("Difference in function value: %g\n", std::abs(g_p_i - g_pp_i)); 
        if (std::abs(g_p_i - g_pp_i)>=thresh ) {
          printf("Function value after %i steps %g\n", n_steps, g_pp_i);
          get_oa(i, cs, pp, 0, true);// Generate Tangents at current point
          for (int j = 0; j < n; j++) {
            p[j] = pp[j]; // Move all previous points to the current one
          }
          g_p_i = g_pp_i; 
          nbG[i]++;
        }
      }
      
    }


      
   si.applyCuts(cs);
   delete [] p;
   delete [] pp;
   delete [] up;
  }
}

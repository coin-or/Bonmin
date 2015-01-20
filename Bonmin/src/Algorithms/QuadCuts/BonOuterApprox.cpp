// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 10/16/2007
#include "BonOuterApprox.hpp"
#include "BonBabSetupBase.hpp"
#include "BonOsiTMINLPInterface.hpp"
#include "BonTypes.hpp"
#include <vector>
#include <sstream>
namespace Bonmin {

int OuterApprox::nTimesCalled = 0;

void 
OuterApprox::initialize(Bonmin::BabSetupBase & b){
   b.options()->GetNumericValue("tiny_element", tiny_, "bonmin.");
   b.options()->GetNumericValue("very_tiny_element", veryTiny_, "bonmin.");
}
void
OuterApprox::extractLinearRelaxation(Bonmin::OsiTMINLPInterface &minlp,
                        OsiSolverInterface *si, 
                        const double * x, bool getObj)
  {
  
    int n;
    int m;
    int nnz_jac_g;
    int nnz_h_lag;
    Ipopt::TNLP::IndexStyleEnum index_style;

    Bonmin::TMINLP2TNLP * model = minlp.problem();
    
    //Get problem information
    model->get_nlp_info( n, m, nnz_jac_g, nnz_h_lag, index_style);
  
    vector<int> iRow(nnz_jac_g);
    vector<int> jCol(nnz_jac_g);
    vector<double> vals(nnz_jac_g);
  
    //get Jacobian
    model->eval_jac_g(n, x, 1, m, nnz_jac_g, &iRow[0], &jCol[0], NULL);
    model->eval_jac_g(n, x, 1, m, nnz_jac_g, NULL, NULL, &vals[0]);
 
    //Put jacobian arrays in c style 
    if(index_style == Ipopt::TNLP::FORTRAN_STYLE){
      for(int k = 0 ; k < nnz_jac_g ; k++){
       iRow[k]--;
       jCol[k]--;
      }
    }
    vector<double> g(m);
    model->eval_g(n, x, 1, m, &g[0]);
  
    vector<double> rowLow(m);
    vector<double> rowUp(m);
    vector<double> colUp(n);
    vector<double> colLow(n);

    model->get_bounds_info(n, &colLow[0], &colUp[0], m, &rowLow[0], &rowUp[0]);

    double infty = si->getInfinity();
    double nlp_infty = infty;
    
   vector<Ipopt::TNLP::LinearityType> const_types(m); 
   model->get_constraints_linearity(m, &const_types[0]);
    for(int i = 0 ; i < m ; i++) {
       if(rowLow[i] > - nlp_infty){
         rowLow[i] = (rowLow[i] - g[i]);
         if(1 || const_types[i] != Ipopt::TNLP::LINEAR) rowLow[i] -= 1e-07;
       }
        if(rowUp[i] < nlp_infty){
          rowUp[i] =  (rowUp[i] - g[i]);
         if(1 || const_types[i] != Ipopt::TNLP::LINEAR) rowUp[i] += 1e-07;
        }
    }
    
    //Then convert everything to a CoinPackedMatrix
    //Go through values, clean coefficients and fix bounds
    for(int i = 0 ; i < nnz_jac_g ; i++) {
      if(const_types[iRow[i]] != Ipopt::TNLP::LINEAR || //For linear just copy is fine.
         cleanNnz(vals[i],colLow[jCol[i]], colUp[jCol[i]],
                  rowLow[iRow[i]], rowUp[iRow[i]],
                  x[jCol[i]],
                  rowLow[iRow[i]],
                  rowUp[iRow[i]], tiny_, veryTiny_)) {      
            rowLow[iRow[i]] += vals[i] * x[jCol[i]];
            rowUp[iRow[i]] += vals[i] *x[jCol[i]];
         }
    }
    CoinPackedMatrix mat(true, &iRow[0], &jCol[0], &vals[0], nnz_jac_g);
    mat.setDimensions(m,n); // In case matrix was empty, this should be enough
    
    vector<double> obj(n,0.);
    
    si->loadProblem(mat, &colLow[0], &colUp[0], &obj[0], &rowLow[0], &rowUp[0]);

    for(int i = 0 ; i < n ; i++) {
      if(minlp.isInteger(i))
        si->setInteger(i);
    }
    if(getObj) {
       if(model->hasLinearObjective()){
         std::cout<<"Linear stuff"<<std::endl;
         double zero;
         model->eval_f(n, &obj[0], 1, zero);
         si->setDblParam(OsiObjOffset, -zero);
         //if(fabs(zero - 0) > 1e-10)
           //addObjVar = true;
         //else { 
           //Copy the linear objective and don't create a dummy variable.
           model->eval_grad_f(n, x, 1,&obj[0]);
           si->setObjective(&obj[0]);
         //}
      }
     else {
        //add variable alpha
        //(alpha should be empty in the matrix with a coefficient of -1 and unbounded)
        CoinPackedVector a;
        si->addCol(a,-si->getInfinity(), si->getInfinity(), 1.);
    
        // Now get the objective cuts
        // get the gradient, pack it and add the cut
        model->eval_grad_f(n, x, 1,&obj[0]);
        double ub;
        model->eval_f(n, x, 1, ub);
        ub*=-1;
        double lb = -1e300;
        CoinPackedVector objCut;
        CoinPackedVector * v = &objCut;
        v->reserve(n);
        for(int i = 0; i<n ; i++) {
         if(nnz_jac_g)
         {
          if(cleanNnz(obj[i],colLow[i], colUp[i],
              -si->getInfinity(), 0,
              x[i],
              lb,
              ub, tiny_, veryTiny_)) {
            v->insert(i,obj[i]);
            lb += obj[i] * x[i];
            ub += obj[i] * x[i];
          }
         }
         else //Unconstrained problem can not put clean coefficient
         {
             if(cleanNnz(obj[i],colLow[i], colUp[i],
              -si->getInfinity(), 0,
              x[i],
              lb,
              ub, 1e-03, 1e-08)) {
            v->insert(i,obj[i]);
            lb += obj[i] * x[i];
            ub += obj[i] * x[i];
             }
         }
        }
      v->insert(n,-1);
      si->addRow(objCut, lb, ub);
      }
    }

#if 0
    std::ostringstream os;
    os<<"OA_"<<nTimesCalled;
    nTimesCalled++;
    std::string f_name = os.str();
    si->writeMps(f_name.c_str());
#endif
  }
  
}


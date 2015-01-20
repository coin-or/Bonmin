// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
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


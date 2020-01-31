// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 10/06/2007

#include "BonTMINLP2Quad.hpp"
#include <climits>

using namespace Ipopt;

//#define DEBUG
namespace Bonmin {

    TMINLP2TNLPQuadCuts::TMINLP2TNLPQuadCuts(const SmartPtr<Bonmin::TMINLP> tminlp):
      TMINLP2TNLP(tminlp)
     {
       // Fill the locally stored hessian matrix
       
       // Get the number of nonzoeroes in the matrix
       const int nnz_h = TMINLP2TNLP::nnz_h_lag();
       curr_nnz_jac_ = TMINLP2TNLP::nnz_jac_g();
       if(nnz_h > 0){
         int * jCol = new int [nnz_h];
         int * iRow = new int [nnz_h];
         
         TMINLP2TNLP::eval_h(num_variables(), NULL, false, 
                             0., TMINLP2TNLP::num_constraints(), NULL, false, 
                             nnz_h, jCol, iRow, NULL);
  
         for(int i = 0 ; i < nnz_h ; i++){
#ifndef NDEBUG
           bool inserted = 
#endif
                    H_.insert(std::make_pair( std::make_pair(jCol[i], iRow[i]), 
                              std::make_pair(i, -1))).second;
           assert(inserted == true);
         }
         delete [] jCol;
         delete [] iRow;
       }
       assert(nnz_h == (int) H_.size());
       obj_.reserve(TMINLP2TNLP::num_variables());
     }


    /** Copy Constructor 
      * \warning source and copy point to the same tminlp_.
      */
    TMINLP2TNLPQuadCuts::TMINLP2TNLPQuadCuts(const TMINLP2TNLPQuadCuts &other):
      TMINLP2TNLP(other),
      quadRows_(other.quadRows_),
      H_(),
      curr_nnz_jac_(other.curr_nnz_jac_),
      obj_(other.obj_)
      {
       // Get the number of nonzoeroes in the matrix
       const size_t nnz_h = TMINLP2TNLP::nnz_h_lag();

       if(nnz_h > 0){
         int * jCol = new int [nnz_h];
         int * iRow = new int [nnz_h];
         int m = TMINLP2TNLP::num_constraints() - (int)quadRows_.size(); 
         TMINLP2TNLP::eval_h(num_variables(), NULL, false, 
                             0., m, NULL, false, 
                             (int)nnz_h, jCol, iRow, NULL);

         for(size_t i = 0 ; i < nnz_h ; i++){
#ifndef NDEBUG
           bool inserted = 
#endif
                    H_.insert(std::make_pair( std::make_pair(jCol[i], iRow[i]), 
                              std::make_pair((int)i, -1))).second;
           assert(inserted == true);
         }
         delete [] jCol;
         delete [] iRow;
        }
         assert(nnz_h == H_.size());

        //Properly create quadRows_
       for(size_t i = 0 ; i < quadRows_.size() ; i++){
         quadRows_[i] = new QuadRow(*quadRows_[i]);
        }

	int offset = TMINLP2TNLP::index_style() == Ipopt::TNLP::FORTRAN_STYLE;
        for(unsigned int i = 0 ; i < quadRows_.size() ; i++){
         quadRows_[i]->add_to_hessian(H_, offset);
        }
      }

    
    /** Destructor */
    TMINLP2TNLPQuadCuts::~TMINLP2TNLPQuadCuts(){
      for(unsigned int i = 0 ; i < quadRows_.size() ; i++){
         delete quadRows_[i];
      }
    }

    
     bool TMINLP2TNLPQuadCuts::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
        Index& nnz_h_lag,
        TNLP::IndexStyleEnum& index_style){
        bool ret_val = TMINLP2TNLP::get_nlp_info(n,m,nnz_jac_g, nnz_h_lag, index_style);
        nnz_h_lag = (int)H_.size();
        nnz_jac_g = curr_nnz_jac_;
        //printf("Dinmension in TMINLP2Quad are %i\n", curr_nnz_jac_);
        return ret_val;
      }

    /** This call is just passed onto parent class and add bounds of quadratic
        cuts*/
     bool TMINLP2TNLPQuadCuts::get_bounds_info(Index n, Number* x_l, Number* x_u,
        Index m, Number* g_l, Number* g_u){
        bool ret_val = TMINLP2TNLP::get_bounds_info(n, x_l, x_u, 
                     m, g_l, g_u);
        return ret_val;
      }

    bool 
    TMINLP2TNLPQuadCuts::get_constraints_linearity(Index m, LinearityType* const_types)
    {
      bool ret_val = TMINLP2TNLP::get_constraints_linearity(m - (int)quadRows_.size(), const_types);
      const_types += m - (int)quadRows_.size();
      for(unsigned int i = 0 ; i < quadRows_.size() ; i++){
        if(quadRows_[i]->isLinear())
          const_types[i] = TNLP::LINEAR;
        else
          const_types[i] = TNLP::NON_LINEAR;
      }
      return ret_val;
    }

    /** This call is just passed onto parent class and add 
        lambda for quadratic cuts*/
     bool TMINLP2TNLPQuadCuts::get_starting_point(Index n, bool init_x, Number* x,
        bool init_z, Number* z_L, Number* z_U,
        Index m, bool init_lambda,
        Number* lambda){
         return TMINLP2TNLP::get_starting_point(n, init_x, x, init_z, z_L, z_U, m, init_lambda, lambda);
    }

    /** Method that returns scaling parameters (passed to parent all quadratic
        not scaled). 
     */
     bool TMINLP2TNLPQuadCuts::get_scaling_parameters(Number& obj_scaling,
                                        bool& use_x_scaling, Index n,
                                        Number* x_scaling,
                                        bool& use_g_scaling, Index m,
                                        Number* g_scaling){
           assert(num_constraints() == m);
           bool retval = get_scaling_parameters(obj_scaling, use_x_scaling, n, x_scaling, use_g_scaling, m - (int)quadRows_.size(), g_scaling);
           if(use_g_scaling){
             g_scaling += m - (int)quadRows_.size();
             CoinFillN(g_scaling, (int)quadRows_.size(), 1.);}
           return retval;
      }

  /** Returns the value of the objective function in x*/
  bool 
  TMINLP2TNLPQuadCuts::eval_f(Index n, const Number* x, bool new_x,
        Number& obj_value){
    if(obj_.empty()){
       return TMINLP2TNLP::eval_f(n, x, new_x, obj_value);
    }
    if(new_x){
       TMINLP2TNLP::eval_f(n,x, new_x, obj_value);
    }
    obj_value = c_;
    assert(n == (int) obj_.size());
    for(int i = 0 ; i < n ; i++){
      obj_value += obj_[i] * x[i];
    }
    return true;
  }

  /** Returns the vector of the gradient of
    *  the objective w.r.t. x */
  bool 
  TMINLP2TNLPQuadCuts::eval_grad_f(Index n, const Number* x, bool new_x,
        Number* grad_f){
    if(obj_.empty()){
      return TMINLP2TNLP::eval_grad_f(n, x, new_x, grad_f);}
    if(new_x){
      TMINLP2TNLP::eval_grad_f(n, x, new_x, grad_f);}
    assert(n == (int) obj_.size());
    for(int i = 0 ; i < n ; i++){
      grad_f[i] = obj_[i];
    }
   return true;
  }

  bool TMINLP2TNLPQuadCuts::eval_gi(Index n, const Number* x, bool new_x,
                           Index i, Number& gi)
  {
    int m_orig = num_constraints() - (int)quadRows_.size();
    if(i < m_orig){
       return TMINLP2TNLP::eval_gi(n, x, new_x, i, gi);
    }
    i -= m_orig;
     gi = quadRows_[i]->eval_f(x, new_x);
     return false;
  }


    /** Returns the vector of constraint values in x (appends constraint values for quadratics).*/
     bool TMINLP2TNLPQuadCuts::eval_g(Index n, const Number* x, bool new_x,
        Index m, Number* g){
       int m_tminlp = m - (int)quadRows_.size();
       bool retval = TMINLP2TNLP::eval_g(n, x, new_x, m_tminlp, g);
       g+= (m_tminlp);
       for(unsigned int i = 0 ; i < quadRows_.size() ; i++){
         g[i] = quadRows_[i]->eval_f(x, new_x);
       }
      return retval;
    }

    /** Returns the jacobian of the
     *  constraints. The vectors iRow and jCol only need to be set
     *  once. The first call is used to set the structure only (iRow
     *  and jCol will be non-NULL, and values will be NULL) For
     *  subsequent calls, iRow and jCol will be NULL. */
     bool TMINLP2TNLPQuadCuts::eval_jac_g(Index n, const Number* x, bool new_x,
        Index m, Index nele_jac, Index* iRow,
        Index *jCol, Number* values){
        int n_ele_orig =  TMINLP2TNLP::nnz_jac_g();
        int m_orig = m - (int)quadRows_.size();
	int offset = TMINLP2TNLP::index_style() == Ipopt::TNLP::FORTRAN_STYLE;

        bool retval = TMINLP2TNLP::eval_jac_g(n, x, new_x, m_orig ,
                                n_ele_orig, iRow, jCol, values);
        if(values == NULL){
           assert(iRow != NULL);
           assert(jCol != NULL);
           iRow += n_ele_orig;
           jCol += n_ele_orig;
           for(unsigned int i = 0 ; i < quadRows_.size() ; i++){
             const int & nnz = quadRows_[i]->nnz_grad();
             Ipopt::Index mi = m_orig + i + offset;
             CoinFillN(iRow, nnz, mi);
             quadRows_[i]->gradiant_struct(nnz, jCol, offset);
             iRow += nnz;
             jCol += nnz;
           }
         }
         else {
           assert(iRow == NULL);
           assert(jCol == NULL);
           values += n_ele_orig;
           for(unsigned int i = 0 ; i < quadRows_.size() ; i++){
             const int & nnz = quadRows_[i]->nnz_grad();
             quadRows_[i]->eval_grad(nnz, x, new_x, values);
             values+=nnz;
            }
          }
       return retval;
    }

  bool TMINLP2TNLPQuadCuts::eval_grad_gi(Index n, const Number* x, bool new_x,
                                Index i, Index& nele_grad_gi, Index* jCol,
                                Number* values)
  {
    int m_orig = num_constraints() - (int)quadRows_.size();
    if(i < m_orig){
       return TMINLP2TNLP::eval_grad_gi(n, x, new_x, i, nele_grad_gi, jCol, values);
    }
    i -= m_orig;
    int offset = TMINLP2TNLP::index_style() == Ipopt::TNLP::FORTRAN_STYLE;
    if(values == NULL){
      assert(jCol != NULL);
      nele_grad_gi = quadRows_[i]->nnz_grad();
      quadRows_[i]->gradiant_struct(nele_grad_gi, jCol, offset);
    }
    else{
      assert(jCol == NULL);
      quadRows_[i]->eval_grad(nele_grad_gi, x, new_x, values);
    }
    return false;
  }
    /** Return the hessian of the
     *  lagrangian. The vectors iRow and jCol only need to be set once
     *  (during the first call). The first call is used to set the
     *  structure only (iRow and jCol will be non-NULL, and values
     *  will be NULL) For subsequent calls, iRow and jCol will be
     *  NULL. This matrix is symmetric - specify the lower diagonal
     *  only */
     bool TMINLP2TNLPQuadCuts::eval_h(Index n, const Number* x, bool new_x,
        Number obj_factor, Index m, const Number* lambda,
        bool new_lambda, Index nele_hess,
        Index* iRow, Index* jCol, Number* values){
        if(!obj_.empty()) obj_factor = 0;
        if(values == NULL){
           assert(iRow != NULL);
           assert(jCol != NULL);
#ifdef DEBUG
           std::cout<<"Hessian structure"<<std::endl;
#endif
	   int nnz = 0;
           int nnz_h_lag_orig = TMINLP2TNLP::nnz_h_lag();
           int nnz_sup = nnz_h_lag_orig;
           for(AdjustableMat::iterator i = H_.begin() ; i != H_.end() ; i++){
              if(i->second.second == -1){
                 assert(i->second.first < nnz_h_lag_orig);
               }
               else {
                 assert(i->second.second > 0);
                 assert(i->second.first >= nnz_h_lag_orig);
                 i->second.first = nnz_sup;
                 nnz_sup++;
               }
              iRow[i->second.first] = i->first.first;
              jCol[i->second.first] = i->first.second;
#ifdef DEBUG
              printf("iRow %i, jCol %i : nnz %i\n",
                     i->first.second, i->first.first, 
                     i->second.first);
#endif
              //assert(*jCol >= *iRow);
              nnz++;
           }
	   assert(nnz == (int) H_.size());
           return true;
         }
         else {
#ifdef DEBUG
           std::cout<<"Computing hessian"<<std::endl;
#endif
           assert(iRow == NULL);
           assert(jCol == NULL);
           int nnz_h_lag_orig = TMINLP2TNLP::nnz_h_lag();
           int m_orig = m - (int)quadRows_.size();
           bool ret_val = TMINLP2TNLP::eval_h(n, x, new_x, obj_factor, m_orig, lambda, new_lambda,
                            nnz_h_lag_orig, iRow, jCol, values);
	   CoinZeroN(values + nnz_h_lag_orig, (int)H_.size() - nnz_h_lag_orig);
           for(unsigned int i = 0 ; i < quadRows_.size() ; i++){
             quadRows_[i]->eval_hessian(lambda[i + m_orig], values);
            }
            return ret_val;
          }
      }
    //@}

  /** Method to add quadratic cuts .*/
  void 
  TMINLP2TNLPQuadCuts::addCuts(const Cuts & cuts, bool safe){
     assert(cuts.sizeColCuts() == 0);
#ifdef DEBUG
     printf("Adding %i cuts\n", cuts.sizeRowCuts());
#endif
     int offset = TMINLP2TNLP::index_style() == Ipopt::TNLP::FORTRAN_STYLE;
     
     g_l_.reserve(g_l_.size() + cuts.sizeQuadCuts() + cuts.sizeRowCuts());
     g_u_.reserve(g_u_.size() + cuts.sizeQuadCuts() + cuts.sizeRowCuts());
     quadRows_.reserve(quadRows_.size() + cuts.sizeQuadCuts() + cuts.sizeRowCuts());

     int n = cuts.sizeQuadCuts();
     for(int i = 0 ; i < n ; i++){
       g_l_.push_back(cuts.quadCut(i).lb());
       g_u_.push_back(cuts.quadCut(i).ub());
       quadRows_.push_back(new QuadRow(cuts.quadCut(i)));
       quadRows_.back()->add_to_hessian(H_, offset);
       curr_nnz_jac_ += quadRows_.back()->nnz_grad();
     }
     addRowCuts((OsiCuts) cuts, safe); 
     duals_sol_.resize(g_l_.size() + 2*x_l_.size(), 0.);
     x_init_.resize(g_l_.size() + 3*x_l_.size(), 0.);
     duals_init_ = x_init_() + x_l_.size();
  }

  /** Method to add array of OsiRowCut figuring out which is quadratic (slow!)..*/
  void TMINLP2TNLPQuadCuts::addCuts(unsigned int numcuts, 
                                    const OsiRowCut ** cuts){
#ifdef DEBUG
     printf("Adding %i cuts\n", numcuts);
#endif
     int offset = TMINLP2TNLP::index_style() == Ipopt::TNLP::FORTRAN_STYLE;
     g_l_.reserve(g_l_.size() + numcuts);
     g_u_.reserve(g_u_.size() + numcuts);
     quadRows_.reserve(quadRows_.size() + numcuts);
     for(unsigned int i = 0 ; i < numcuts ; i++){
       g_l_.push_back(cuts[i]->lb());
       g_u_.push_back(cuts[i]->ub());

       const QuadCut * quadCut = dynamic_cast<const QuadCut *> (cuts[i]); 
       if(quadCut){
         quadRows_.push_back(new QuadRow(*quadCut));
         quadRows_.back()->add_to_hessian(H_, offset);
       }
       else 
        quadRows_.push_back(new QuadRow(*cuts[i]));
       curr_nnz_jac_ += quadRows_.back()->nnz_grad();
     }
     duals_sol_.resize(g_l_.size() + 2*x_l_.size(), 0.);
     x_init_.resize(g_l_.size() + 3*x_l_.size(), 0.);
     duals_init_ = x_init_() + x_l_.size();
  } 
  /** Method to add OsiCuts figuring out which is quadratic (slow!)..*/
  void TMINLP2TNLPQuadCuts::addCuts(const OsiCuts& cuts){
     assert(cuts.sizeColCuts() == 0);
#ifdef DEBUG
     printf("Adding %i cuts\n", cuts.sizeRowCuts());
#endif

     const Cuts * quadCuts = dynamic_cast<const Cuts *>(&cuts);
     if(quadCuts) {
        addCuts(*quadCuts, true);
        return;}

     addRowCuts(cuts, true);
  } 
  /** Method to add OsiCuts eventualy figuring out which is quadratic.*/
  void TMINLP2TNLPQuadCuts::addRowCuts(const OsiCuts& cuts, bool safe){
    // Check with rowCuts are quadratics and move them to quadratic cuts.
    int n = cuts.sizeRowCuts(); 
     g_l_.reserve(g_l_.size() + n);
     g_u_.reserve(g_u_.size() + n);
     quadRows_.reserve(quadRows_.size() + n);

    int offset = TMINLP2TNLP::index_style() == Ipopt::TNLP::FORTRAN_STYLE;

    for(int i = 0 ; i < n ; i++){
      g_l_.push_back(cuts.rowCut(i).lb());
      g_u_.push_back(cuts.rowCut(i).ub());
      if(safe == false){
      assert(dynamic_cast<const QuadCut *> (cuts.rowCutPtr(i)) == NULL);
      }
      else {
       const QuadCut * cut = dynamic_cast<const QuadCut *> (cuts.rowCutPtr(i)); 
       if(cut){
         quadRows_.push_back(new QuadRow(*cut));
         quadRows_.back()->add_to_hessian(H_, offset);
         curr_nnz_jac_ += quadRows_.back()->nnz_grad();
         continue;
      } 
    } 
    quadRows_.push_back(new QuadRow(cuts.rowCut(i)));
    curr_nnz_jac_ += quadRows_.back()->nnz_grad();
  }
     duals_sol_.resize(g_l_.size() + 2*x_l_.size(), 0.);
     x_init_.resize(g_l_.size() + 3*x_l_.size(), 0.);
     duals_init_ = x_init_() + x_l_.size();
 }

  /** Method which removes a set of cuts.*/
  void TMINLP2TNLPQuadCuts::removeCuts(unsigned int n,const int * idxs){
     if(n == 0) return;
     vector< int > order(quadRows_.size());
     int m_tminlp = num_constraints() - (int)quadRows_.size();
      //delete the pointers
       for(unsigned int k = 0; k < n ; k++){//Erase
       int idx = idxs[k] - m_tminlp ;
       quadRows_[idx]->remove_from_hessian(H_);
       curr_nnz_jac_ -= quadRows_[idx]->nnz_grad();
       delete quadRows_[idx];
       quadRows_[idx] = NULL;}

     for(unsigned int i = 0 ; i < order.size() ; i++){
        order[i] = i;
     }
     for(unsigned int i = 0 ; i < n ; i++){
        assert(idxs[i] - m_tminlp >= 0);
        order[ idxs[i] - m_tminlp ] = INT_MAX;
    } 

    std::sort(order.begin(), order.end());


    int i;
    double * g_l = g_l_() + m_tminlp;
    double * g_u = g_u_() + m_tminlp;
    for(i = 0 ; order[i] < INT_MAX ; i++){
      assert(order[i] >= i);
      quadRows_[i] = quadRows_[order[i]]; 
      g_l[i] = g_l[order[i]];
      g_u[i] = g_u[order[i]];
    }
    quadRows_.erase(quadRows_.begin() + i, quadRows_.end());
    g_l_.erase(g_l_.begin() + m_tminlp + i, g_l_.end());
    g_u_.erase(g_u_.begin() + m_tminlp + i, g_u_.end());
 }

void
TMINLP2TNLPQuadCuts::printH(){
  int nnz = 0;
  for(AdjustableMat::iterator i = H_.begin() ; i != H_.end() ; i++){
     std::cout<<"nnz: "<<nnz
	     <<"jCol: "<<i->first.first
	     <<", iRow "<<i->first.second<<std::endl;
    nnz++;
  }
}

void
TMINLP2TNLPQuadCuts::set_linear_objective(int n_var, const double * obj, double c_0){
  assert(n_var == TMINLP2TNLP::num_variables());
  obj_.resize(n_var);
  CoinCopyN(obj, n_var, obj_());
  c_ = c_0;
}
}//Ends Bonmin namespace
 

// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 10/06/2007

#include "BonQuadRow.hpp"
#include <cfloat>
//#define DEBUG
namespace Bonmin{

QuadRow::QuadRow():
  c_(0),
  a_(),
  Q_(),
  grad_evaled_(false)
{
}    

QuadRow::QuadRow(const QuadRow &other):
  c_(other.c_),
  a_(other.a_),
  Q_(other.Q_),
  g_(),
  a_grad_idx_(),
  Q_row_grad_idx_(),
  Q_col_grad_idx_(),
  Q_hessian_idx_(),
  grad_evaled_(false)
{
  initialize();
}

QuadRow & QuadRow::operator=(const QuadRow &rhs){
  if(this != &rhs){
    c_ = rhs.c_;
    a_ = rhs.a_;
    Q_ = rhs.Q_;
    Q_hessian_idx_.clear();
    g_.clear();
    a_grad_idx_.clear();
    Q_row_grad_idx_.clear();
    Q_col_grad_idx_.clear();
    initialize();
    //H_Hes_idx_ = rhs.H_Hes_idx_;
   grad_evaled_ = false;
  }
  return (*this);
}

QuadRow::QuadRow(const QuadCut &cut):
  c_(0),
  a_(cut.row()),
  Q_(cut.Q(), cut.type())
  {
    initialize(); 
  }

QuadRow& 
QuadRow::operator=(const QuadCut &cut){
    c_ = cut.c();
    a_ = cut.row();
    Q_ = cut.Q();
    Q_.make_upper_triangular(cut.type());
    g_.clear();
    a_grad_idx_.clear();
    Q_row_grad_idx_.clear();
    Q_col_grad_idx_.clear();
    //Q_hessian_idx.empty();
    //H_Hes_idx_.empty()
    initialize();
    return (*this);
}


QuadRow::QuadRow(const OsiRowCut &cut):
  c_(0),
  a_(cut.row()),
  Q_()
  {
    initialize(); 
  }

QuadRow& 
QuadRow::operator=(const OsiRowCut &cut){
    c_ = 0;
    a_ = cut.row();
    Q_ = TMat();
    g_.empty();
    a_grad_idx_.empty();
    Q_row_grad_idx_.clear();
    Q_col_grad_idx_.clear();
    //Q_hessian_idx.empty();
    //H_Hes_idx_.empty()
    initialize();
    return (*this);
}

void
QuadRow::initialize(){
    //Check that Q_ is upper triangular
    for(int i = 0 ; i < Q_.nnz_ ; i++){
      assert(Q_.jCol_[i] >= Q_.iRow_[i]);}
    grad_evaled_ = false;

   int n;
    // Construct a map to store the non-zero elements of the gradient.
   n = a_.getNumElements();
   a_grad_idx_.reserve(n);

   // Put the linear elements
   const int * indices = a_.getIndices();
   const double * elems = a_.getElements();

   for(int i = 0 ; i < n ; i++){
    std::pair<gStore::iterator, bool> res = g_.insert(std::make_pair(indices[i], std::make_pair(elems[i],0.)));
    a_grad_idx_.push_back(res.first);
    }
   // Put the quadratics first rows
   n = Q_.numNonEmptyRows();
   const TMat::RowS& nonEmptyRows = Q_.nonEmptyRows();
   Q_row_grad_idx_.reserve(n);

   for(TMat::RowS::const_iterator i = nonEmptyRows.begin() ; i != nonEmptyRows.end() ; i++){
    std::pair<gStore::iterator, bool> res = g_.insert(std::make_pair(i->first, std::make_pair(0.,0.)));
    Q_row_grad_idx_.push_back(res.first);
   }

   //Now columns
   n = Q_.numNonEmptyCols();
   const TMat::RowS& nonEmptyCols = Q_.nonEmptyCols();
   Q_col_grad_idx_.reserve(n);

   for(TMat::RowS::const_iterator i = nonEmptyCols.begin() ; i != nonEmptyCols.end() ; i++){
    std::pair<gStore::iterator, bool> res = g_.insert(std::make_pair(i->first, std::make_pair(0.,0.)));
    Q_col_grad_idx_.push_back(res.first);
   }

}

/** Print quadratic constraint.*/
void
QuadRow::print(){
  std::cout<<"constant term "<<c_<<std::endl;
  int * a_ind = a_.getIndices();
  const double * a_el = a_.getElements();
  int n = a_.getNumElements();
  std::cout<<"Linear term (size "<<n<<"): ";
  for(int i = 0 ; i < n ; i++)
  {
    std::cout<<a_el[i]<<" * x["<<a_ind[i]<<"]\t";
    if(i && ( (i % 5) == 0 ) ) std::cout<<std::endl<<"\t\t";
  }
}
/** Evaluate quadratic form.*/
double 
QuadRow::eval_f(const double *x, bool new_x){
  //if(new_x){
    internal_eval_grad(x);//}
  double value = c_;// Constant

  //Linear part
  int * a_ind = a_.getIndices();
  const double * a_el = a_.getElements();
  int n = a_.getNumElements();
  for(int i = 0 ; i < n ; i++)
  {
    value += a_el[i] * x[a_ind[i]];
    
  }

  //Quadratic part
  for(gStore::iterator i = g_.begin() ; i != g_.end() ; i++){
    value += i->second.second * x[i->first];
  }
  return value;
}

/** Get number of non-zeroes in the gradiant.*/
int 
QuadRow::nnz_grad(){
  return static_cast<int>(g_.size());}
/** Get structure of gradiant */
void 
QuadRow::gradiant_struct(const int nnz, int * indices, bool offset){
  int n = 0;
  for(gStore::iterator i = g_.begin() ; i != g_.end() ; i++){
    indices[n++] = i->first + offset;
  }
  assert(n == nnz);
  assert(nnz == (int) g_.size());
}

/** Evaluate gradiant of quadratic form.*/
void 
QuadRow::eval_grad(const int nnz, const double * x, bool new_x, double * values){

#ifdef DEBUG
   // Output relevant components of x
   for(gStore::iterator i = g_.begin() ; i != g_.end() ; i++){
     printf("x[%i] = %g,  ",i->first, x[i->first]);
   }
#endif
  //if(new_x){
    internal_eval_grad(x);//}
  int n = 0;
#ifdef DEBUG
  std::cout<<"Computing gradient"<<std::endl;
#endif
  for(gStore::iterator i = g_.begin() ; i != g_.end() ; i++){
#ifdef DEBUG
    printf("%i: %g, %g\n", i->first, i->second.second, i->second.first);
#endif
    values[n++] = 2*i->second.second + i->second.first;
  }
  assert (nnz == (int) g_.size());
}

void
QuadRow::internal_eval_grad(const double *x){
   // Zero out gradiant storage
   for(gStore::iterator i = g_.begin() ; i != g_.end() ; i++){
     i->second.second = 0;
   }


   const TMat::RowS & nonEmptyRows = Q_.nonEmptyRows();
   int k = 0;//Iterates on Q_grad_idx_;
   for(TMat::RowS::const_iterator ii = nonEmptyRows.begin() ; ii != nonEmptyRows.end();
       ii++, k++){
    double value = 0;
    assert(ii->first == Q_.iRow_[Q_.rowOrdering_[ii->second]]);
    for(int i = ii->second ; i < Q_.nnz_ && ii->first == Q_.iRow_[Q_.rowOrdering_[i]] ; i++)
    {
        value += x[Q_.jCol_[Q_.rowOrdering_[i]]] * Q_.value_[Q_.rowOrdering_[i]];
    }
    Q_row_grad_idx_[k]->second.second += value;
    assert(Q_row_grad_idx_[k]->first == ii->first);
   }
   const TMat::RowS & nonEmptyCols = Q_.nonEmptyCols();
   k = 0;//Iterates on Q_grad_idx_;
   for(TMat::RowS::const_iterator ii = nonEmptyCols.begin() ; ii != nonEmptyCols.end();
       ii++, k++){
    double value = 0;
    assert(ii->first == Q_.jCol_[Q_.columnOrdering_[ii->second]]);
    for(int i = ii->second ; i < Q_.nnz_ && ii->first == Q_.jCol_[Q_.columnOrdering_[i]] ; i++)
    {
      if(Q_.iRow_[Q_.columnOrdering_[i]] != Q_.jCol_[Q_.columnOrdering_[i]])
        value += x[Q_.iRow_[Q_.columnOrdering_[i]]] * Q_.value_[Q_.columnOrdering_[i]];
    }
    Q_col_grad_idx_[k]->second.second += value;
    assert(Q_col_grad_idx_[k]->first == ii->first);
   }

   grad_evaled_ = true;
}

void
QuadRow::add_to_hessian(AdjustableMat &H, bool offset){
  assert(Q_hessian_idx_.empty());
  for(int i = 0 ; i < Q_.nnz_ ; i++){
     std::pair<int, int> e;
     e = std::make_pair(Q_.jCol_[i] + offset, Q_.iRow_[i] + offset);
     AdjustableMat::iterator pos = H.find(e);
     if(pos != H.end()){//Already exists
       if(pos->second.second != -1)
          pos->second.second++;
        Q_hessian_idx_.push_back(pos); 
     }
     else {
        std::pair<AdjustableMat::iterator, bool> res = 
            H.insert(std::make_pair(e, std::make_pair((int)H.size(), 1)));
        assert(res.second == true);
        Q_hessian_idx_.push_back(res.first);
     }
  } 
}

void
QuadRow::remove_from_hessian(AdjustableMat &H){
  for(int i = 0 ; i < Q_.nnz_ ; i++){
     if(Q_hessian_idx_[i]->second.second != -1)
        Q_hessian_idx_[i]->second.second--;
     if(Q_hessian_idx_[i]->second.second == 0){
        H.erase(Q_hessian_idx_[i]);
     }
  }
  Q_hessian_idx_.clear();
}

/** Return hessian values (i.e. Q_) in values.*/
 void 
 QuadRow::eval_hessian(double lambda, double * values){
  for(int i = 0 ; i < Q_.nnz_ ; i++){
#ifdef DEBUG
     printf("iRow %i, jCol %i, value %g , nnz %i\n",
            Q_hessian_idx_[i]->first.second,
            Q_hessian_idx_[i]->first.first,
            Q_.value_[i],
            Q_hessian_idx_[i]->second.first);
#endif
     values[Q_hessian_idx_[i]->second.first] += (lambda * 2 * Q_.value_[i]);
  }
}

}// Ends Bonmin namespace


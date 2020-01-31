// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 10/06/2007

#include "BonTMatrix.hpp"

namespace Bonmin{

/** Copy constructor.*/
TMat::TMat(const TMat &other):
  iRow_(NULL), jCol_(NULL), value_(NULL), nnz_(other.nnz_),
  capacity_(other.nnz_), columnOrdering_(other.columnOrdering_),
  rowOrdering_(other.rowOrdering_), nonEmptyRows_(), nonEmptyCols_(){
   iRow_ = CoinCopyOfArray(other.iRow_, other.nnz_);
   jCol_ = CoinCopyOfArray(other.jCol_, other.nnz_);
   value_ = CoinCopyOfArray(other.value_, other.nnz_);
  }

/** Construct from a CoinPackedMatrix*/
TMat::TMat(const CoinPackedMatrix &M, MatrixStorageType T):
  iRow_(NULL), jCol_(NULL), value_(NULL), nnz_(M.getNumElements()),
  capacity_(M.getNumElements()), columnOrdering_(), rowOrdering_(),
  nonEmptyRows_(), nonEmptyCols_(){
  create(M);
  make_upper_triangular(T);
}
/** Assignment operator.*/
TMat& 
TMat::operator=(const TMat &rhs){
  if(this != &rhs){
    freeSpace();
    nnz_ = rhs.nnz_;
    capacity_ = rhs.capacity_;
    iRow_ = CoinCopyOfArray(rhs.iRow_, rhs.nnz_);
    jCol_ = CoinCopyOfArray(rhs.jCol_, rhs.nnz_);
    value_ = CoinCopyOfArray(rhs.value_, rhs.nnz_);
    columnOrdering_ = rhs.columnOrdering_;
    rowOrdering_ = rhs.rowOrdering_; 
    nonEmptyCols_.clear();
    nonEmptyRows_.clear();
  }
  return (*this);
}

/** Assignment from a CoinPackedMatrix.*/
TMat & 
TMat::operator=(const CoinPackedMatrix &M){
  freeSpace();
  columnOrdering_.clear();
  rowOrdering_.clear();
  nnz_ = capacity_ = M.getNumElements();
  create(M); 
  return (*this);
}

void TMat::create(const CoinPackedMatrix &M){ 
  // Allocate arrays;
  iRow_ = new int[capacity_];
  jCol_ = new int[capacity_];
  value_ = new double[capacity_];

  int * iRow = iRow_;
  int * jCol = jCol_;
  if(!M.isColOrdered()){// Have to swap
    std::cout<<"Matrix is not col ordered"<<std::endl;
    iRow = jCol_;
    jCol = iRow_;
  }
  
  // Now we can safely assume that M is colorderd.
  int numcols = M.getMajorDim();
  const int * start = M.getVectorStarts();
  const int * length = M.getVectorLengths();
  const int * indice = M.getIndices();
  const double * value = M.getElements();
  int nnz = 0;
  for(int i = 0 ; i < numcols ; i++){
    int begin = start[i];
    int end = start[i] + length[i];
    for(int k = begin ; k < end ; k++){
      value_[nnz] = value[k];
      iRow[nnz] = indice[k];
      jCol[nnz++] = i;
    }
  }
  assert(nnz==nnz_);
}

// Destructor
TMat::~TMat(){
   delete [] iRow_;
   delete [] jCol_;
   delete [] value_;
}

/** Put the non-empty rows of quadratic form currently stored into
    indices. Allocate the array and returns its size.*/
int 
TMat::numNonEmptyRows(){
  if(nnz_ == 0) return 0;
  orderByRows();
  nonEmptyRows_.clear();
  nonEmptyRows_.push_back(std::pair<int, int>(iRow_[rowOrdering_[0]], 0));
  int num = 1;
  for(int i = 1 ; i < nnz_ ; i++){
    if(iRow_[rowOrdering_[i]] > nonEmptyRows_.back().first){
      nonEmptyRows_.push_back(std::pair<int, int>(iRow_[rowOrdering_[i]],i));
      num++;
    }
  }
  return num;
}


/** Put the non-empty rows of quadratic form currently stored into
    indices. Allocate the array and returns its size.*/
int 
TMat::numNonEmptyCols(){
  if(nnz_ == 0) return 0;
  orderByColumns();
  nonEmptyCols_.clear();
  nonEmptyCols_.push_back(std::pair<int, int>(jCol_[columnOrdering_[0]], 0));
  int num = 1;
  for(int i = 1 ; i < nnz_ ; i++){
    if(jCol_[columnOrdering_[i]] > nonEmptyCols_.back().first){
      nonEmptyCols_.push_back(std::pair<int, int>(jCol_[columnOrdering_[i]],i));
      num++;
    }
  }
  return num;
}
/** Remove the duplicated entries.*/
void 
TMat::removeDuplicates(){
  orderByRows();
  int j = 0;
  for(int i = 1; i < nnz_ ; i++){
    if((jCol_[rowOrdering_[i]] == jCol_[rowOrdering_[j]]) && 
       (iRow_[rowOrdering_[i]] == iRow_[rowOrdering_[j]])){
       value_[rowOrdering_[j]] += value_[rowOrdering_[i]];
    }
    else{
      jCol_[rowOrdering_[++j]] = jCol_[rowOrdering_[i]];
      iRow_[rowOrdering_[j]] = iRow_[rowOrdering_[i]];
      value_[rowOrdering_[j]] = value_[rowOrdering_[i]];
    }
  }
  resizeAndCopyArray(jCol_, j, capacity_);
  resizeAndCopyArray(iRow_, j, capacity_);
  resizeAndCopyArray(value_, j, capacity_);
  nnz_ = j;
}

void
TMat::make_upper_triangular(const MatrixStorageType &T){
   switch (T){
     case Upper:
       for(int i = 0 ; i < nnz_ ; i++){
          assert(jCol_[i] >= iRow_[i]);
       }
       break;
     case Lower:
       for(int i = 0 ; i < nnz_ ; i++){
          assert(jCol_[i] <= iRow_[i]);
       }
       make_lower_to_be_upper();
       break;
     case Full:
       make_full_upper_triangular();
       break;
   }
   for(int i = 0 ; i < nnz_ ; i++){
      assert(jCol_[i] >= iRow_[i]);
   }
}

/** Assuing that this is representing the lower triangle of a symetric matrix
   makes it the upper triangle.*/
void
TMat::make_lower_to_be_upper(){
  int * buff = iRow_;
  iRow_ = jCol_;
  jCol_ = buff;
}

/** Assuming that this is representing a quadratic form. Makes it upper diagonal.*/
void
TMat::make_full_upper_triangular(){
  // Make it upper triangular
  for(int i = 0 , j = 0; i < nnz_ && j < nnz_ ;){
    if(iRow_[i] < jCol_[i]){//swap the two entries
      int buf = iRow_[i];
      iRow_[i] = jCol_[i];
      jCol_[i] = buf;
    }
  }
  // add up the duplicated entries
  removeDuplicates();

  //Now divide all non-diagonal entries by two;
  for(int i = 0 ; i < nnz_ ; i++){
    if(jCol_[i] == iRow_[i]){//skip
      continue;
    }
    assert(iRow_[i] < jCol_[i]);
    value_[i] /= 2.;
  }
}

}//Ends Bonmin namepace


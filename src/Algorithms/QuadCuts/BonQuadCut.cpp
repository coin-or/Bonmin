// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 10/06/2007

#include "BonQuadCut.hpp"

namespace Bonmin {

  QuadCut::QuadCut():
    OsiRowCut(),
    c_(0),
    Q_(),
    type_(Upper)
    {}

  QuadCut::QuadCut(const QuadCut & other):
    OsiRowCut(other),
    c_(other.c_),
    Q_(other.Q_),
    type_(other.type_){
    }

  QuadCut &
  QuadCut::operator=(const QuadCut & rhs){
    if(this != &rhs){
       OsiRowCut::operator=(rhs);
       c_ = rhs.c_;
       Q_ = rhs.Q_;
       type_ = rhs.type_;
    }
   return *this;
  }

  OsiRowCut *
  QuadCut::clone() const {
    return new QuadCut(*this);
   }

  QuadCut::~QuadCut(){
  }

  /// Compute cut violation
  double 
  QuadCut::violated(const double * solution) const{
    double rhs = c_;
    rhs += row().dotProduct(solution);
    const int * indice = Q_.getIndices();
    const double * val = Q_.getElements();
    const int * start = Q_.getVectorStarts();
    const int * length = Q_.getVectorLengths();
    int n = Q_.getMajorDim();

    for(int i = 0 ; i < n ; i++){
    	int s=start[i];
	int l=length[i];
      for(int k = s ; k < s+l ; k++){
        if(i == indice[k]) rhs += solution[i] * solution[indice[k]] * val[k];
        else rhs += 2*solution[i] * solution[indice[k]] * val[k];
      }
    }

#if 0
    for(int i = 0 ; i < n ; i++){
      for(int k = 0 ; k < length[i] ; k++){
        if(i == indice[k]) rhs += solution[i] * solution[indice[k]] * val[k];
        else rhs += 2*solution[i] * solution[indice[k]] * val[k];
      }
      start += length[i];
    }
#endif
    return std::max(lb() - rhs, rhs - ub());
  }
    


  void
  QuadCut::print() const{
    std::cout<<"Quadratic cut has lower bound "<<lb()<<" and upper bound "<<ub()
             <<std::endl;
    
    std::cout<<"Linear part has "<<row().getNumElements()<<" non zeroes:"
             <<std::endl;  

    const int& nElem = row().getNumElements();
    const int * indices = row().getIndices();
    const double * elements = row().getElements();

    for(int i = 0 ; i < nElem ; i++){
      if(i > 0 && elements[i] > 0.)
         std::cout<<"+ ";
        std::cout<< elements[i] <<" x["<<indices[i]<<"]\t";
       if(i > 0 && i % 5 == 0) std::cout<<std::endl;
    }
    std::cout<<std::endl;
    if(Q_.getNumElements()){
      std::cout<<"Quadratic part is given by the matrix:"<<std::endl;
      Q_.dumpMatrix();
    }
   }
  Cuts::Cuts():
    OsiCuts(),
    quadCuts_(0){
   } 

  Cuts::Cuts(const Cuts & other):
    OsiCuts(other),
    quadCuts_(other.quadCuts_.size()){
    for(unsigned int i = 0 ; i < quadCuts_.size() ; i++){
      quadCuts_[i] = new QuadCut(*other.quadCuts_[i]);
    }
  }

  Cuts &
  Cuts::operator=(const Cuts & rhs){
    if(this != &rhs){
      OsiCuts::operator=(rhs);
      for(unsigned int i = 0 ; i < quadCuts_.size() ; i++)
      {
         delete quadCuts_[i];
      }
      quadCuts_.resize(rhs.quadCuts_.size());
      for(unsigned int i = 0 ; i < quadCuts_.size() ; i++){
        quadCuts_[i] = new QuadCut(*rhs.quadCuts_[i]);
      }
    }
    return *this;
  }

  Cuts::~Cuts(){
          for(unsigned int i = 0 ; i < quadCuts_.size() ; i++)
      {
         delete quadCuts_[i];
      }
   }
      
void
Cuts::printCuts() const {
  OsiCuts::printCuts();
  std::cout<<quadCuts_.size()<<" quadratic cuts."<<std::endl;
  for(unsigned int i = 0 ; i < quadCuts_.size() ; i++){
    quadCuts_[i]->print();
  }
}

} /* Ends Bonmin namespace.*/

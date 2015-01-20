// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 10/16/2007
#ifndef SepaSepaTMINLP2OsiLP_H
#define SepaSepaTMINLP2OsiLP_H

#include <cmath>
#include <cstdio>
#include "IpSmartPtr.hpp"
#include "IpTNLP.hpp"
#include "BonTypes.hpp"
#include "BonTMINLP2OsiLP.hpp"


namespace Sepa {

  /** A transformer class to build outer approximations i.e. transfomrs nonlinear programs into linear programs.*/
  class SepaTMINLP2OsiLP: public Bonmin::TMINLP2OsiLP {

  public:

   /** Default constructor.*/
   SepaTMINLP2OsiLP():
    TMINLP2OsiLP(),
    num_approx_(0)
   {}

   /** Copy constructor.*/
   SepaTMINLP2OsiLP(const SepaTMINLP2OsiLP & other):
    TMINLP2OsiLP(other), 
    num_approx_(other.num_approx_){
    }

   TMINLP2OsiLP * clone() const{
     return new SepaTMINLP2OsiLP(*this);
   }

   /** Assignment operator.*/
   SepaTMINLP2OsiLP & operator=(const SepaTMINLP2OsiLP& rhs){
    if(this != & rhs){
      TMINLP2OsiLP::operator=(rhs);
      num_approx_ = rhs.num_approx_;
    }
    return (*this);
   }

   void set_num_approx(int v){
     num_approx_ = v;
   }

   /** Destructor.*/
   ~SepaTMINLP2OsiLP(){}

   /** Build the Outer approximation of model_ in x and put it in si.*/
   virtual void extract(OsiSolverInterface *si, 
                const double * x, bool getObj) ;

   /** Get OAs of nonlinear constraints in x.*/
   virtual void get_oas(OsiCuts & cs, 
                const double * x, bool getObj, bool global) const;

   virtual void get_refined_oa(OsiCuts & cs) const;
   /** Get OA of one constraints in x.*/
   virtual void get_oa(int iRow, OsiCuts & cs, 
                const double * x, bool getObj, bool global) const;

   void add_outer_description(OsiSolverInterface &si) ;
   void add_outer_description_function_values(OsiSolverInterface &si) ;
   private:
   int num_approx_;
  };


}

#endif


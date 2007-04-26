// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/18/2007
#ifndef BonCouenneSetup_H
#define BonCouenneSetup_H
#include "BonBabSetupBase.hpp"

namespace Bonmin{
  class CouenneSetup : public BabSetupBase{
public:
    /** Default constructor*/
    CouenneSetup():
    BabSetupBase(){}
    
    /** Copy constructor.*/
    CouenneSetup(const CouenneSetup& other):
      BabSetupBase(other){}
    
    /** virtual copy constructor.*/
    virtual BabSetupBase * clone() const{
      return new CouenneSetup(*this);
    }

    /** Initialize from command line arguments.*/
    void InitializeBonmin(char **& argv);
    /** register the options */
    virtual void registerOptions();
    /** Register all Couenne options.*/
    static void registerAllOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions);
    
    /** Get the basic options if don't already have them.*/
    virtual void readOptionsFile(){
      if(readOptions_) return;
      BabSetupBase::readOptionsFile("couenne.opt");}
  };
  
}
#endif

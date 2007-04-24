// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/18/2007
#include "BonAmplSetup.hpp"

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
    
    /** Register all Couenne options.*/
    static void registerAllOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions);
    
    /** Get the basic options if don't already have them.*/
    virtual void defaultBasicOptions(){
      if(GetRawPtr(options_) != NULL && GetRawPtr(roptions_) != NULL &&  GetRawPtr(journalist_) != NULL) return;
      BasicSetup b;
      Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions = b.roptions();
      CouenneSetup::registerAllOptions(roptions);
      b.Initialize("couenne.opt");
      options_ = b.options();
      roptions_ = b.roptions();
      journalist_ = b.journalist();}
  };
  
}


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
//AW #include "BonBabSetupBase.hpp"
#include "BonBonminSetup.hpp"

struct ASL;



namespace Bonmin{
  
  class SmartAsl : public Ipopt::ReferencedObject{
public:
    ASL * asl;
    SmartAsl():
      Ipopt::ReferencedObject(),
      asl(NULL)
    {}
    virtual ~SmartAsl();
  };
  
  class CouenneSetup : public BonminSetup{
public:
    /** Default constructor*/
    CouenneSetup():
    BonminSetup(),
    aslfg_(NULL){}
    
    /** Copy constructor.*/
    CouenneSetup(const CouenneSetup& other):
      BonminSetup(other),
      aslfg_(NULL){}
    
    /** virtual copy constructor.*/
    virtual BabSetupBase * clone() const{
      return new CouenneSetup(*this);
    }
    
    virtual ~CouenneSetup();
    /** Initialize from command line arguments.*/
    void InitializeCouenne(char **& argv);
    /** register the options */
    virtual void registerOptions();
    /** Register all Couenne options.*/
    static void registerAllOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions);
    
    /** Get the basic options if don't already have them.*/
    virtual void readOptionsFile(){
      if(readOptions_) return;
      BabSetupBase::readOptionsFile("couenne.opt");}

    /** Get the cutoff value from the initial solve */
    double getCutOff() const {
      return cutoff_;
    }
private:
      SmartPtr<SmartAsl> aslfg_;

    // Cutoff value after initialSolve
    double cutoff_;
  };
  
}
#endif

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
#include "CbcFeasibilityBase.hpp"

struct ASL;

class CouenneCutGenerator;


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
  
  class CouenneSetup : public BabSetupBase{
public:
    /** Default constructor*/
    CouenneSetup():
    BabSetupBase(),
    aslfg_(NULL),
    CouennePtr_ (NULL),
    displayStats_ (false) {}

    /** Copy constructor.*/
    CouenneSetup(const CouenneSetup& other):
      BabSetupBase(other),
      aslfg_(NULL),
      displayStats_ (other.displayStats_) {}
    
    /** virtual copy constructor.*/
    virtual BabSetupBase * clone() const{
      return new CouenneSetup(*this);
    }
    
    /// destructor
    virtual ~CouenneSetup();

    /** Initialize from command line arguments. */
    void InitializeCouenne(char **& argv);

    /** register the options */
    virtual void registerOptions();
    /** Register all Couenne options.*/
    static void registerAllOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions);
    
    /** Get the basic options if don't already have them.*/
    virtual void readOptionsFile(){
      if (readOptions_) return;
      BabSetupBase::readOptionsFile("couenne.opt");
    }

    /// return pointer to cut generator (used to get pointer to problem)
    CouenneCutGenerator *couennePtr () const
    {return CouennePtr_;}

    /// true if one wants to display statistics at the end of program
    bool displayStats () 
    {return displayStats_;}

    /// add cut generators
    void addMilpCutGenerators ();

    /// modify parameter (used for MaxTime)
    inline void setDoubleParameter (const DoubleParameter &p, const double val)
    {doubleParam_ [p] = val;}

    /// modify parameter (used for MaxTime)
    inline double getDoubleParameter (const DoubleParameter &p) const
    {return doubleParam_ [p];}

private:
      SmartPtr<SmartAsl> aslfg_;

    /// hold a reference to Couenne cut generator to delete it at
    /// last. The alternative would be to clone it every time the
    /// CouenneSolverInterface is cloned (that is, at each call of
    /// Optimality-based bound tightening).
    CouenneCutGenerator *CouennePtr_;

    /// true if one wants to display statistics at the end of program
    bool displayStats_;
  };
}

#endif

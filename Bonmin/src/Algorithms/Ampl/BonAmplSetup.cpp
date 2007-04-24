// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/15/2007

#include "BonAmplSetup.hpp"
namespace Bonmin{
  void BonminAmplSetup::initializeBonmin(char **& argv){
    /* Get the basic options. */
    BasicSetup b;
    Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions = b.roptions();
    BonminSetup::registerAllOptions(roptions);
    b.Initialize("bonmin.opt");

    /* Read the model.*/
    SmartPtr<AmplTMINLP> model = new AmplTMINLP(ConstPtr(b.journalist()), b.options(), argv, NULL, "bonmin", NULL);
    b.mayPrintDoc();
    setBasicOptions(b);
    BonminSetup::initializeBonmin(GetRawPtr(model));
  }
  
   void BonminAmplSetup::initializeBonmin(AmplInterface &toFill, char **& argv){
    /* Get the basic options. */
    BasicSetup b;
    Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions = b.roptions();
    BonminSetup::registerAllOptions(roptions);
    b.Initialize("bonmin.opt");
    
    /* Read the model.*/
    SmartPtr<AmplTMINLP> model = new AmplTMINLP(ConstPtr(b.journalist()), b.options(), argv, NULL, "bonmin", NULL);
    b.mayPrintDoc();
    toFill.initialize(b, GetRawPtr(model));
    
    BonminSetup::initializeBonmin(toFill);
  }
  
  /** Usefull for Bcp */
   void BonminAmplSetup::fillOsiInterface(AmplInterface &toFill, char ** &argv, std::string & options,
                                                           std::string & nl){
   BasicSetup b;
   Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions = b.roptions();
   BonminSetup::registerAllOptions(roptions);
   b.InitializeFromLongString(options);
   
   /* Read the model.*/
   SmartPtr<AmplTMINLP> model = new AmplTMINLP(ConstPtr(b.journalist()), 
                                               b.options(), 
                                               argv, NULL, "bonmin", 
                                               &nl);
   b.mayPrintDoc();   
   toFill.initialize(b, GetRawPtr(model));
}
  
}


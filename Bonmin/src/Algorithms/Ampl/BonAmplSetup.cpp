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
  void BonminAmplSetup::initialize(char **& argv){
    readOptionsFile();
    /* Read the model.*/
    SmartPtr<AmplTMINLP> model = new AmplTMINLP(ConstPtr(journalist()), options(), argv, NULL, "bonmin", NULL);
    mayPrintDoc();
    BonminSetup::initialize(GetRawPtr(model), true);
  }
  
   void 
   BonminAmplSetup::initialize(AmplInterface &toFill, char **& argv){
    Ipopt::SmartPtr<TNLPSolver> solver = toFill.solver();
    setOptionsAndJournalist(solver->roptions(), 
                            solver->options(),
                            solver->journalist());
    /* Get the basic options. */
     readOptionsFile(); 
    /* Read the model.*/
    SmartPtr<AmplTMINLP> model = new AmplTMINLP(ConstPtr(journalist()), options(), argv, NULL, "bonmin", NULL);
    mayPrintDoc();
    toFill.initialize(roptions_, options_, journalist_, GetRawPtr(model));
    BonminSetup::initialize(toFill, true);
  }
 
  /** initialize bonmin with ampl model using the command line arguments reading options and nl file from strings.*/ 
  void 
  BonminAmplSetup::initialize(char **& argv, std::string& opt_file_content, std::string& nl_file_content, bool createContinuousSolver /*= false*/){
    /* Get the basic options. */
    readOptionsString(opt_file_content);
    /* read nl file by creating AmplTMINLP.*/
    SmartPtr<AmplTMINLP> model = new AmplTMINLP(ConstPtr(journalist()), options(), argv, NULL, "bonmin", &nl_file_content);
    mayPrintDoc();
    BonminSetup::initialize(GetRawPtr(model), createContinuousSolver);}
  
  
  /** initialize bonmin with ampl model using the command line arguments and an existing OsiTMINLPInterface reading options and nl file from strings.*/
  void 
  BonminAmplSetup::initialize(AmplInterface &toFill, char **& argv, std::string& opt_file_content, 
                                    std::string& nl_file_content, bool createContinuousSolver /*=  false*/
  ){
    /* Get the basic options. */
    readOptionsString(opt_file_content);
    /* read nl file by creating AmplTMINLP.*/
    SmartPtr<AmplTMINLP> model = new AmplTMINLP(ConstPtr(journalist()), options(), argv, NULL, "bonmin", &nl_file_content);
    mayPrintDoc();
    toFill.initialize(roptions_, options_, journalist_, GetRawPtr(model));
    BonminSetup::initialize(toFill, createContinuousSolver);    
  }
  
  /** Usefull for Bcp */
   void BonminAmplSetup::fillOsiInterface(AmplInterface &toFill, char ** &argv, std::string & options,
                                                           std::string & nl, bool createContinuousSolver /*=  false*/){

     /* Get the basic options. */
     readOptionsString(options);  
   /* Read the model.*/
   SmartPtr<AmplTMINLP> model = 
     new AmplTMINLP(ConstPtr(journalist_), 
                    options_, 
                    argv, NULL, "bonmin", &nl);
   toFill.initialize(roptions(), options_, journalist(), GetRawPtr(model));
}
  
}


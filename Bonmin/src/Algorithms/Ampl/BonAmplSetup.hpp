// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/15/2007

#ifndef BonAmplSetup_H
#define BonAmplSetup_H
#include "BonBonminSetup.hpp"
#include "BonAmplInterface.hpp"

namespace Bonmin
{
  class BonminAmplSetup: public BonminSetup
  {
  public:
    /** initialize bonmin with ampl model using the command line arguments.*/
    void initialize(char **& argv);
    /** initialize bonmin with ampl model using the command line arguments and an existing OsiTMINLPInterface.*/
    void initialize(AmplInterface &toFill, char **& argv);
    /** initialize bonmin with ampl model using the command line arguments reading options and nl file from strings.*/
    void initialize(char **& argv, std::string& opt_file_content, std::string& nl_file_content, bool createContinuousSolver /*= false*/);
    /** initialize bonmin with ampl model using the command line arguments and an existing OsiTMINLPInterface reading options and nl file from strings.*/
    void initialize(AmplInterface &toFill, char **& argv, std::string& opt_file_content, std::string& nl_file_content, bool createContinuousSolver = true);
    /** For Bcp. Initialize the passed OsiTMINLP interface with ampl model using the options and nl files contained in two strings.*/
    void fillOsiInterface(AmplInterface &toFill, char **& argv, std::string & options, std::string & nl, bool createContinuousSolver = true);
  };
}
#endif

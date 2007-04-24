// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/15/2007

#include "BonBonminSetup.hpp"
#include "BonAmplInterface.hpp"

namespace Bonmin{
  class BonminAmplSetup: public BonminSetup {
public:
    /** initialize bonmin with ampl model using the command line arguments.*/ 
    void initializeBonmin(char **& argv);
    /** initialize bonmin with ampl model using the command line arguments and an existing OsiTMINLPInterface.*/
    void initializeBonmin(AmplInterface &toFill, char **& argv);
    /** For Bcp. Initialize the passed OsiTMINLP interface with ampl model using the options and nl files contained in two strings.*/
    static 
      void fillOsiInterface(AmplInterface &toFill, char **& argv, std::string & options, std::string & nl);
};
}


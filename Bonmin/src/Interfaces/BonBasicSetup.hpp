// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/13/2007

#ifndef BonBasicSetup_H
#define BonBasicSetup_H


#include "IpRegOptions.hpp"
#include "IpOptionsList.hpp"

#include "BonBaseOptions.hpp"

#include <string>
#include <fstream>
#include <sstream>


namespace Bonmin{
  /** A class to get the options at the beginning of initialization.*/
  class BasicSetup{
public:
    /** Default constructor.*/
    BasicSetup();
    /** Register the two options.*/
    static void registerOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions>);
    /** Get solver choosen in options.*/
    BaseOptions::Solver getSolver();
    /** Get solver choosen in options.*/
    BaseOptions::Algorithm getAlgorithm();
   /** Initialize options with option file.*/
    void Initialize(const std::string& optFile);
    /** Initialize options with stream.*/
    void Initialize(std::istream& is);
    /** Initialize options with string.*/
    void InitializeFromLongString(std::string & string);
    /** May print the documentation if user asked.*/
    void mayPrintDoc();
    /** Acces storage of Journalist for output */
    Ipopt::SmartPtr<Ipopt::Journalist> journalist(){ return journalist_;}
    
    /** Acces list of Options */
    Ipopt::SmartPtr<Ipopt::OptionsList> options(){return options_;}
    
    /** Access registered Options */
    Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions(){return roptions_;}
private:
    /** Storage of Journalist for output */
    Ipopt::SmartPtr<Ipopt::Journalist> journalist_;
    
    /** List of Options */
    Ipopt::SmartPtr<Ipopt::OptionsList> options_;
    
    /** Registered Options */
    Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions_;
  };
}/* end namespace Bonmin.*/

#endif

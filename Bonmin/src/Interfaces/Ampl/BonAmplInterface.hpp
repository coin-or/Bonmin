// (C) Copyright International Business Machines Corporation and
// Carnegie Mellon University 2004, 2007
//
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
// Andreas Waechter, International Business Machines Corporation
//
// Date : 12/01/2004

#ifndef BonminAmplInterface_H
#define BonminAmplInterface_H
#include "BonOsiTMINLPInterface.hpp"
#include "BonAmplTMINLP.hpp"

class BM_lp;
namespace Bonmin
{
  /** Class for providing an Osi interface to Ipopt with an ampl nl file as input. */
  class AmplInterface: public OsiTMINLPInterface
  {
  public:
    /** Default constructor */
    /** Default constructor only available for Bonmin's friends and child classes.*/
    AmplInterface();
    /**@name Methods to input a problem */
    //@{
    virtual void readAmplNlFile(char **& argv, Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
        Ipopt::SmartPtr<Ipopt::OptionsList> options,
        Ipopt::SmartPtr<Ipopt::Journalist> journalist,
        std::string* nl_file_content  = NULL);
    //@}
    /** Copy constructor */
    AmplInterface(const AmplInterface &other);
    /// Clone
    virtual OsiSolverInterface * clone(bool CopyData = true);

    /// Destructor
    virtual ~AmplInterface();


    /** Fast access to AmplTMINLP */
    const AmplTMINLP * amplModel() const
    {
      return GetRawPtr(amplTminlp_);
    }
    /** To set some application specific defaults. */
    virtual void setAppDefaultOptions(Ipopt::SmartPtr<Ipopt::OptionsList> Options);

  protected:
    /** Read variables and row names in .col and .row files.*/
    void readNames() ;

    /** TMINLP problem (the original problem usually an AmplTMINLP).*/
    Ipopt::SmartPtr<Bonmin::AmplTMINLP> amplTminlp_;

  private:
    /** Write the ampl solution file or write a bonmin one?*/
    int writeAmplSolFile_;
  };
}
#endif

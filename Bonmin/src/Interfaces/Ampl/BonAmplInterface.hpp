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

namespace Bonmin
{
  /** Class for providing an Osi interface to Ipopt with an ampl nl file as input. */
  class AmplInterface: public OsiTMINLPInterface
  {
  public:
    /** Default constructor */
    AmplInterface();
    /** Constructor with inputed ampl command line (reads model from nl file). */
    AmplInterface(char **& amplArgs);
    /** Constructor with inputed ampl command line (reads model from nl file). */
    AmplInterface(char **& amplArgs, Ipopt::SmartPtr<TNLPSolver> app);
    /** Copy constructor */
    AmplInterface(const AmplInterface &other);
    /// Clone
    virtual OsiSolverInterface * clone(bool CopyData = true);

    /// Destructor
    virtual ~AmplInterface();

    /**@name Methods to input a problem */
    //@{
    /** Read an ampl . nl file from the given filename */
    virtual void readAmplNlFile(char**& filename,
        std::string* ipopt_file_content = NULL,
        std::string* nl_file_content = NULL
        );
    /** write ampl solution file */
    void writeAmplSolFile(std::string message,const double * primalSol = NULL,const double * dualSol = NULL);
    //@}

    /** Fast access to AmplTMINLP */
    const AmplTMINLP * amplModel() const
    {
      return GetRawPtr(amplTminlp_);
    }
    /** To set some application specific defaults. */
    virtual void setAppDefaultOptions(Ipopt::SmartPtr<Ipopt::OptionsList> Options);

  protected:

    /** TMINLP problem (the original problem usually an AmplTMINLP).*/
    Ipopt::SmartPtr<Bonmin::AmplTMINLP> amplTminlp_;
  };
}
#endif

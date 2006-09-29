// (C) Copyright International Business Machines Corporation, Carnegie Mellon University 2004
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
// Carl D. Laird, Carnegie Mellon University,
// Andreas Waechter, International Business Machines Corporation
//
// Date : 12/01/2004


#ifndef IpoptInterface_H
#define IpoptInterface_H


#include "BonOsiTMINLPInterface.hpp"
#include "BonIpoptSolver.hpp"

namespace Bonmin{
/**
   This is class provides an Osi interface for Ipopt
   (so that we can use it for example as the continuous solver in Cbc).
*/

class IpoptInterface : public OsiTMINLPInterface
{
  friend class BonminParam;

public:



class UnsolvedIpoptError: public OsiTMINLPInterface::UnsolvedError
{
 public:
  UnsolvedIpoptError(int errorNum = 10000):
  OsiTMINLPInterface::UnsolvedError(errorNum)
  {}
  virtual const std::string& errorName() const;
 
  virtual const std::string& solverName() const; 
  virtual ~UnsolvedIpoptError(){}
  private:
  static std::string errorNames [17];
  static std::string solverName_;
};

  virtual UnsolvedError * newUnsolvedError(int num){
    return new UnsolvedIpoptError(num);}
  //#############################################################################


  /**@name Constructors and destructors */
  //@{
  /// Default Constructor
  IpoptInterface();
  /** Constructor with given (user) TMINLP.
    \warning In this constructor option file is not read, use readOptionFile to read one.
  */
  IpoptInterface (Ipopt::SmartPtr<Bonmin::TMINLP> tminlp);

  /// Clone
  virtual OsiSolverInterface * clone(bool CopyData=true) const;

  /** Copy constructor
  */
  IpoptInterface (const IpoptInterface &);

  /// Assignment operator
  IpoptInterface & operator=(const IpoptInterface& rhs);

  /// Destructor
  virtual ~IpoptInterface ();


  void extractInterfaceParams();
  //@}

  //---------------------------------------------------------------------------

  /// Return status of last optimization
  Ipopt::ApplicationReturnStatus getOptStatus() const;




  //---------------------------------------------------------------------------
  /**@name WarmStart related methods (those should really do nothing for the moment)*/
  //@{

  /*! \brief Get an empty warm start object

  This routine returns an empty CoinWarmStartBasis object. Its purpose is
  to provide a way to give a client a warm start basis object of the
  appropriate type, which can resized and modified as desired.
  */
  CoinWarmStart *getEmptyWarmStart () const;

  /** Get warmstarting information */
  virtual CoinWarmStart* getWarmStart() const;

  /** Set warmstarting information. Return true/false depending on whether
      the warmstart information was accepted or not. */
  virtual bool setWarmStart(const CoinWarmStart* warmstart);


  void setWarmStartOptions()
  {
    //    app_->Options()->SetIntegerValue("warm_start_init_point", 1);
    app_->Options()->SetStringValue("warm_start_init_point", "yes");
  }
  void unsetWarmStartOptions()
  {
    //app_->Options()->SetIntegerValue("warm_start_init_point", 1);
    app_->Options()->SetStringValue("warm_start_init_point", "no");
    problem_->resetStartingPoint();
  }

  //@}


  //---------------------------------------------------------------------------



  /**@name Control of Ipopt output
   */
  //@{
  void turnOffIpoptOutput();
  void turnOnIpoptOutput();
  //@}

  
//---------------------------------------------------------------------------
protected:

  /** Call Ipopt to solve or resolve the problem and check for errors.*/
  void solveAndCheckErrors(bool doResolve, bool throwOnFailure,
      const char * whereFrom);


  //@}

  /** Warm start strategy :
  <ol>
  <li> no warm start,</li>
  <li> simple warm start (optimal point),</li>
  <li> more elaborate strategies (interior point...).</li>
  </ol>
  */
  int warmStartStrategy_;
  /** flag to say wether options have been printed or not.*/
  static bool hasPrintedOptions;
};
}
#endif

// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2004
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
// Carl D. Laird, Carnegie Mellon University,
// Andreas Waechter, International Business Machines Corporation
//
// Date : 12/01/2004


#include "BonIpoptInterface.hpp"
#include "BonTMINLP.hpp"
#include "BonColReader.hpp"
#include "CoinTime.hpp"
#include "BonIpoptWarmStart.hpp"
#include "BonIpoptInteriorWarmStarter.hpp"
#include <string>
#include <sstream>

#include "IpSolveStatistics.hpp"
#include "BonIpoptSolver.hpp"

using namespace Ipopt;

namespace Bonmin{
bool IpoptInterface::hasPrintedOptions=0;

////////////////////////////////////////////////////////////////////
// Constructors and desctructors                                  //
////////////////////////////////////////////////////////////////////
/// Default Constructor
IpoptInterface::IpoptInterface():
    OsiTMINLPInterface(),
    warmStartStrategy_(1)
{}


/** Constructor with given IpSolver and TMINLP */
IpoptInterface::IpoptInterface (Ipopt::SmartPtr<Bonmin::TMINLP> tminlp):
    OsiTMINLPInterface(),
    warmStartStrategy_(1)
{
  Ipopt::SmartPtr<Bonmin::TNLPSolver> app = new Bonmin::IpoptSolver();
  allocateTMINLP(tminlp, app);
}




/// Clone
OsiSolverInterface * IpoptInterface::clone(bool CopyData) const
{
  //    debugMessage("IpoptInterface::clone(%d)\n", CopyData);
  static int debug_no_clone = 0;
  debug_no_clone++;
  //    assert(debug_no_clone==1);
  return new IpoptInterface(*this);
}

/// Copy constructor
IpoptInterface::IpoptInterface (const IpoptInterface &source):
    OsiTMINLPInterface(source),
    warmStartStrategy_(source.warmStartStrategy_)
{
}




/// Assignment operator
IpoptInterface & IpoptInterface::operator=(const IpoptInterface& rhs)
{
  if(this!= &rhs) {
    OsiTMINLPInterface::operator=(rhs);
  }
  return *this;
}

/// Destructor
IpoptInterface::~IpoptInterface ()
{
}

///////////////////////////////////////////////////////////////////
// WarmStart Information                                                                           //
///////////////////////////////////////////////////////////////////


CoinWarmStart* IpoptInterface::getEmptyWarmStart () const
{
  return (dynamic_cast<CoinWarmStart *>(new IpoptWarmStart(1)));
}


/// Get warmstarting information
CoinWarmStart*
IpoptInterface::getWarmStart() const
{
  if(warmStartStrategy_) {
    if(warmStartStrategy_==2) {
      SmartPtr<IpoptInteriorWarmStarter> warm_starter =
        SmartPtr<IpoptInteriorWarmStarter>(problem_->GetWarmStarter());
      return new IpoptWarmStart(*this, warm_starter);
    }
    else  return new IpoptWarmStart(*this, NULL);
  }
  else
    return new IpoptWarmStart(getNumCols(), getNumRows());
}


bool
IpoptInterface::setWarmStart(const CoinWarmStart* warmstart)
{
  if(!warmstart || !warmStartStrategy_)
    return 0;
  hasBeenOptimized_ = false;
  const IpoptWarmStart * ws = dynamic_cast<const IpoptWarmStart*> (warmstart);
  if(ws->empty())//reset initial point and leave
  {
    unsetWarmStartOptions();
    return 1;
  }
  setWarmStartOptions();
  int numcols = getNumCols();
  int numrows = getNumRows();
  const double * colLow = getColLower();
  const double * colUp = getColUpper();
  for(int i = 0; i < numcols ; i++) {
    CoinWarmStartBasis::Status status = ws->getStructStatus(i);
    if(status == CoinWarmStartBasis::atLowerBound) {
      problem_->setxInit(i,colLow[i]);
      problem_->setDualInit(i + numcols + numrows,0.);
    }
    else if(status == CoinWarmStartBasis::atUpperBound) {
      problem_->setxInit(i,colUp[i]);
      problem_->setDualInit(i + numrows,0.);
    }
    else {
      problem_->setDualInit(i + numrows,0.);
      problem_->setDualInit(i + numcols + numrows, 0.);
    }
  }
  for(int i = 0; i < numrows ; i++) {
    CoinWarmStartBasis::Status status = ws->getArtifStatus(i);
    if(status == CoinWarmStartBasis::atLowerBound) {
      problem_->setDualInit(i,0.);
    }
  }
  int nElem = ws->values()->getNumElements();
  const int * inds = ws->values()->getIndices();
  const double * elems = ws->values()->getElements();

  for(int i = 0 ; i < nElem ; i++) {
    problem_->setxInit(inds[i],elems[i]);
  }

  if(IsValid(ws->warm_starter()))
    problem_->SetWarmStarter(ws->warm_starter());
  return 1;
}


void
IpoptInterface::turnOffIpoptOutput()
{
  std::string opt="print_level";
  app_->Options()->SetIntegerValue(opt, (Index)J_NONE,true);
}
void
IpoptInterface::turnOnIpoptOutput()
{
  std::string opt="print_level";
  app_->Options()->SetIntegerValue(opt, (Index)J_SUMMARY,true);
}



/*******************************************************************************/
// Class for throwing errors reported from Ipopt
/******************************************************************************/

std::string
IpoptInterface::UnsolvedIpoptError::errorNames[17] ={"Solve succeeded",
    "Solved to acceptable level",
    "Infeasible problem detected",
    "Search direction becomes too small",
    "Diverging iterates",
    "User requested stop",
    "Maximum iterations exceeded",
    "Restoration failed",
    "Error in step computation",
    "Not enough degrees of freedom",
    "Invalid problem definition",
    "Invalid option",
    "Invalid number detected",
    "Unrecoverable exception",
    "NonIpopt exception thrown",
    "Insufficient memory",
    "Internal error"};

const std::string &
IpoptInterface::UnsolvedIpoptError::errorName() const
{
  if(errorNum() >=0)
    return errorNames[errorNum()];
  if(errorNum() == -1) return errorNames[6];
  else if(errorNum() == -2) return errorNames[7];
  else if(errorNum() == -3) return errorNames[8];
  else if(errorNum() == -10) return errorNames[9];
  else if(errorNum() == -11) return errorNames[10];
  else if(errorNum() == -12) return errorNames[11];
  else if(errorNum() == -13) return errorNames[12];
  else if(errorNum() == -100) return errorNames[13];
  else if(errorNum() == -101) return errorNames[14];
  else if(errorNum() == -102) return errorNames[15];
  else if(errorNum() == -199) return errorNames[16];
  throw IpoptInterface::SimpleError("UnsolvedError::errorName()","Unrecognized optimization status in ipopt.");
}

std::string IpoptInterface::UnsolvedIpoptError::solverName_ = "Ipopt";

const std::string &
IpoptInterface::UnsolvedIpoptError::solverName() const
{
   return solverName_;
}



/// Return status of last optimization
Ipopt::ApplicationReturnStatus 
IpoptInterface::getOptStatus() const
{
  IpoptSolver* ipoptApp = dynamic_cast<IpoptSolver *>(GetRawPtr(app_));
  return ipoptApp->getOptStatus();
}

void
IpoptInterface::extractInterfaceParams()
{
  OsiTMINLPInterface::extractInterfaceParams();
  app_->Options()->GetEnumValue("warm_start",warmStartStrategy_,"bonmin.");

}

} /* end namespace Bonmin */

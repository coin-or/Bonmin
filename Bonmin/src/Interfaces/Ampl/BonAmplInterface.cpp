// (C) Copyright International Business Machines Corporation and
// Carnegie Mellon University 2004, 2007
//
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
// Andreas Waechter, International Business Machines Corporation
//
// Date : 12/01/2004

#include "BonminConfig.h"

#include "BonAmplInterface.hpp"
#include "BonIpoptSolver.hpp"
#include "BonColReader.hpp"
#ifdef COIN_HAS_FILTERSQP
# include "BonFilterSolver.hpp"
#endif
#include <string>
#include <sstream>

#include "BonTNLP2FPNLP.hpp"

namespace Bonmin
{
  /** Default constructor no initialization*/
  AmplInterface::AmplInterface():
      OsiTMINLPInterface(), amplTminlp_(NULL)
  {}

  /** Copy constructor */
  AmplInterface::AmplInterface(const AmplInterface &other):
      OsiTMINLPInterface(other), amplTminlp_(NULL)
  {
    amplTminlp_ = dynamic_cast<Bonmin::AmplTMINLP *> (GetRawPtr(tminlp_));
  }
/// Clone
  OsiSolverInterface *
  AmplInterface::clone(bool CopyData )
  {
    return new AmplInterface(*this);
  }

///Destructor
  AmplInterface::~AmplInterface()
  {
    amplTminlp_ = NULL;
  }


  void AmplInterface::readAmplNlFile(char **& argv, Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
      Ipopt::SmartPtr<Ipopt::OptionsList> options,
      Ipopt::SmartPtr<Ipopt::Journalist> journalist,
      std::string* nl_file_content /* = NULL*/
                                    )
  {
    if (!IsValid(app_)) {
      createApplication(roptions, options, journalist, "bonmin.");
    }
    // set the default options... expect_infeasible, etc...
    if (!IsValid(tminlp_)) {
      amplTminlp_ = new AmplTMINLP(Ipopt::ConstPtr(app_->journalist()), app_->roptions(), app_->options(), argv,
          NULL, appName() , nl_file_content);
      tminlp_ = GetRawPtr(amplTminlp_);
    }
    else {
      AmplTMINLP * amplTMINLP = dynamic_cast<AmplTMINLP *> (GetRawPtr(tminlp_));
      if (amplTMINLP) {
        AmplTMINLP * newAmpl = amplTMINLP->createEmpty();
        newAmpl->Initialize(ConstPtr(app_->journalist()), app_->roptions(), app_->options(), argv,
            NULL, appName() , nl_file_content);
        amplTminlp_ = newAmpl;
        tminlp_ = GetRawPtr(amplTminlp_);
      }
      else {
        amplTminlp_ = new AmplTMINLP(ConstPtr(app_->journalist()), app_->roptions(), app_->options(), argv,
            NULL, appName() , nl_file_content);
        tminlp_ = GetRawPtr(amplTminlp_);
      }
    }
    problem_ = new TMINLP2TNLP(tminlp_);
    feasibilityProblem_ = new TNLP2FPNLP
        (Ipopt::SmartPtr<Ipopt::TNLP>(Ipopt::GetRawPtr(problem_)));
  if(feasibility_mode_){
    problem_to_optimize_ = GetRawPtr(feasibilityProblem_);
  }
  else {
    problem_to_optimize_ = GetRawPtr(problem_);
  }

    int numcols = getNumCols();
    if (obj_)
      delete [] obj_;
    obj_ = new double[numcols];
    CoinFillN(obj_,numcols,1.);
    setStrParam(OsiProbName, std::string(argv[1]));
    extractInterfaceParams();
    hasBeenOptimized_ = false;
    //Read columns and row names if they exists
    readNames();
  }

  void
  AmplInterface::setAppDefaultOptions(Ipopt::SmartPtr<Ipopt::OptionsList> Options)
{}


  void
  AmplInterface::readNames()
  {
    std::string probName;
    getStrParam(OsiProbName, probName);
    NamesReader colRead(probName, ".col");
    if (colRead.readFile()) {
      OsiNameVec colNames;
      colRead.copyNames(colNames);
      setColNames(colNames, 0, static_cast<int>(colNames.size()), 0);
    }

    NamesReader rowRead(probName, ".row");
    if (rowRead.readFile()) {
      OsiNameVec rowNames;
      rowRead.copyNames(rowNames);
      setRowNames(rowNames, 0, static_cast<int>(rowNames.size()), 0);
    }


  }

} /* end namespace Bonmin. */

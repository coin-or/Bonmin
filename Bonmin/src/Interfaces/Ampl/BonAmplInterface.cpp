#include "BonminConfig.h"

#include "BonAmplInterface.hpp"
#include "BonIpoptSolver.hpp"
#ifdef COIN_HAS_FILTERSQP
# include "BonFilterSolver.hpp"
#endif
#include <string>
#include <sstream>


namespace Bonmin
{
  /** Default constructor */
  AmplInterface::AmplInterface(): OsiTMINLPInterface(), amplTminlp_(NULL)
  {}

  /** Constructor with inputed ampl command line (reads model from nl file)*/
  AmplInterface::AmplInterface(char **& amplArgs,
			       SmartPtr<TNLPSolver> app /* =  NULL */)
      :
      OsiTMINLPInterface(),
      amplTminlp_(NULL)
  {
    readAmplNlFile(amplArgs, app, NULL, NULL);
  }

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

  /** Read an ampl . nl file from the given filename */
  void
  AmplInterface::readAmplNlFile(char**& filename,
      Ipopt::SmartPtr<TNLPSolver> app /* = NULL */,
      std::string* ipopt_file_content /* = NULL */,
      std::string* nl_file_content /* = NULL */)
  {
    // If nothing is given in 'app', we determine from the options
    // which solver is to be used
    if (IsNull(app)) {
      // AW: The following is not a nice solution, since we read
      // everything twice, if FilterSQP was chosen

      //We need to build dummy solver objects to get the options,
      //determine which is the solver to use and register all the
      //options
      app_ = new IpoptSolver();
      OsiTMINLPInterface forOption(GetRawPtr(app_));
      int ival;
      forOption.solver()->Options()->GetEnumValue("nlp_solver", ival,"bonmin.");
      NLPSolverChoice NLPchoice = NLPSolverChoice(ival);

      if(NLPchoice == Bonmin::AmplInterface::FilterSQP) {
#ifdef COIN_HAS_FILTERSQP
	app_ = new Bonmin::FilterSolver();
#else
	std::cerr<<"Bonmin not configured to run with FilterSQP"<<std::endl;
	throw -1;
#endif
      }
      else if (NLPchoice != Bonmin::AmplInterface::Ipopt) {
	std::cerr<<"Trying to use unknown solver."<<std::endl;
	throw -1;
      }
    }
    else {
      app_ = app->clone();
    }

    SmartPtr<RegisteredOptions> roptions = app_->RegOptions();
    register_ALL_options(roptions);

    // Call initalize to open output
    app_->Initialize("");
    // Read the bonmin.opt input file
    if (ipopt_file_content == NULL) {
      app_->Initialize("bonmin.opt");
    }
    else {
      std::stringstream ss(ipopt_file_content->c_str());
      app_->Initialize(ss);
    }


    // set the default options... expect_infeasible, etc...
    if (!IsValid(tminlp_)) {
      amplTminlp_ = new AmplTMINLP(Ipopt::ConstPtr(app_->Jnlst()), app_->Options(), filename,
          NULL, appName() , nl_file_content);
      tminlp_ = GetRawPtr(amplTminlp_);
    }
    else {
      AmplTMINLP * amplTMINLP = dynamic_cast<AmplTMINLP *> (GetRawPtr(tminlp_));
      if (amplTMINLP) {
        AmplTMINLP * newAmpl = amplTMINLP->createEmpty();
        newAmpl->Initialize(ConstPtr(app_->Jnlst()), app_->Options(), filename,
            NULL, appName() , nl_file_content);
        amplTminlp_ = newAmpl;
        tminlp_ = GetRawPtr(amplTminlp_);
      }
      else {
        amplTminlp_ = new AmplTMINLP(ConstPtr(app_->Jnlst()), app_->Options(), filename,
            NULL, appName() , nl_file_content);
        tminlp_ = GetRawPtr(amplTminlp_);
      }
    }
    problem_ = new TMINLP2TNLP(tminlp_);//, *app_->Options());

    bool print_options_documentation;
    app_->Options()->GetBoolValue("print_options_documentation",
        print_options_documentation, "");
    if (print_options_documentation) {
      std::list<std::string> categories;
      categories.push_back("bonmin branch-and-bound options");
      categories.push_back("bonmin options for robustness");
      categories.push_back("bonmin options for non-convex problems");
      categories.push_back("bonmin options : B-Hyb specific options");
//    roptions->OutputLatexOptionDocumentation2(*app_->Jnlst(),categories);
      roptions->OutputOptionDocumentation(*(Ipopt::ConstPtr(app_->Jnlst())),categories);
    }

    int numcols = getNumCols();
    if (obj_)
      delete [] obj_;
    obj_ = new double[numcols];
    CoinFillN(obj_,numcols,1.);
    setStrParam(OsiProbName, std::string(filename[1]));
    extractInterfaceParams();
    hasBeenOptimized_ = false;
    feasibilityProblem_ = new TNLP2FPNLP
        (Ipopt::SmartPtr<TNLP>(Ipopt::GetRawPtr(problem_)));
  }

  /** write ampl solution file */
  void
  AmplInterface::writeAmplSolFile(std::string message,const double * primalSol,const double * dualSol)
  {
    TMINLP * tminlp = GetRawPtr(tminlp_);
    AmplTMINLP * ampl_tminlp = dynamic_cast<AmplTMINLP *> (tminlp);
    if (ampl_tminlp)
      ampl_tminlp->write_solution(message,primalSol,dualSol);
    else
      std::cerr<<"Errot can not write .sol file for non ampl problem"<<std::endl;
  }

void 
AmplInterface::setAppDefaultOptions(Ipopt::SmartPtr<Ipopt::OptionsList> Options)
{}


} /* end namespace Bonmin. */

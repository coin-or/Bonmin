#include "BonAmplInterface.hpp"
#include <string>
#include <sstream>


namespace Bonmin
{
  /** Default constructor */
  AmplInterface::AmplInterface(): IpoptInterface(), amplTminlp_(NULL)
  {}

  /** Constructor with inputed ampl command line (reads model from nl file)*/
  AmplInterface::AmplInterface(char **& amplArgs)
      :
      IpoptInterface(),
      amplTminlp_(NULL)
  {
    readAmplNlFile(amplArgs, NULL, NULL);
  }

  /** Copy constructor */
  AmplInterface::AmplInterface(const AmplInterface &other):
      IpoptInterface(other), amplTminlp_(NULL)
  {
    amplTminlp_ = dynamic_cast<Bonmin::AmplTMINLP *> (GetRawPtr(tminlp_));
  }
/// Clone
  AmplInterface *
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
      std::string* ipopt_file_content,
      std::string* nl_file_content)
  {



    app_ = new IpoptSolver;

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
    IpoptSolver * ipopt = dynamic_cast<IpoptSolver *> (GetRawPtr(app_));
    if (!IsValid(tminlp_)) {
      amplTminlp_ = new AmplTMINLP(ConstPtr(ipopt->getIpoptApp().Jnlst()), app_->Options(), filename,
          NULL, appName() , nl_file_content);
      tminlp_ = GetRawPtr(amplTminlp_);
    }
    else {
      AmplTMINLP * amplTMINLP = dynamic_cast<AmplTMINLP *> (GetRawPtr(tminlp_));
      if (amplTMINLP) {
        AmplTMINLP * newAmpl = amplTMINLP->createEmpty();
        newAmpl->Initialize(ConstPtr(ipopt->getIpoptApp().Jnlst()), app_->Options(), filename,
            NULL, appName() , nl_file_content);
        amplTminlp_ = newAmpl;
        tminlp_ = GetRawPtr(amplTminlp_);
      }
      else {
        amplTminlp_ = new AmplTMINLP(ConstPtr(ipopt->getIpoptApp().Jnlst()), app_->Options(), filename,
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
      roptions->OutputOptionDocumentation(*(ipopt->getIpoptApp().Jnlst()),categories);
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

}

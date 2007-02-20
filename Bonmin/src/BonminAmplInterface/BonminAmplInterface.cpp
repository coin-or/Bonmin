#include "BonminAmplInterface.hpp"
#include <string>
#include <sstream>

/** Default constructor */
BonminAmplInterface::BonminAmplInterface(): IpoptInterface(), amplTminlp_(NULL)
{}
  
/** Constructor with inputed ampl command line (reads model from nl file)*/ 
BonminAmplInterface::BonminAmplInterface(char **& amplArgs, bool create_stdout)
:
IpoptInterface(),
amplTminlp_(NULL)
{
  readAmplNlFile(amplArgs, NULL, NULL, create_stdout);
}

/** Copy constructor */
BonminAmplInterface::BonminAmplInterface(const BonminAmplInterface &other):
          IpoptInterface(other), amplTminlp_(NULL)
{
  amplTminlp_ = dynamic_cast<Ipopt::AmplTMINLP *> (GetRawPtr(tminlp_));
}
/// Clone
BonminAmplInterface * 
BonminAmplInterface::clone(bool CopyData )
{
  return new BonminAmplInterface(*this); 
}

///Destructor
BonminAmplInterface::~BonminAmplInterface()
{amplTminlp_ = NULL;}

/** Read an ampl . nl file from the given filename */
void
BonminAmplInterface::readAmplNlFile(char**& filename,
    std::string* ipopt_file_content,
    std::string* nl_file_content,
    bool create_stdout)
{



  app_ = new Ipopt::IpoptApplication( create_stdout);

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
  set_ipopt_minlp_default(app_->Options());

  if(!IsValid(tminlp_)) {
        amplTminlp_ = new AmplTMINLP(ConstPtr(app_->Jnlst()), app_->Options(), filename,
        NULL, appName() , nl_file_content);
        tminlp_ = GetRawPtr(amplTminlp_);
  }
  else {
    AmplTMINLP * amplTMINLP = dynamic_cast<AmplTMINLP *> (GetRawPtr(tminlp_));
    if(amplTMINLP) {
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
  problem_ = new Ipopt::TMINLP2TNLP(tminlp_, *app_->Options());

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
    roptions->OutputOptionDocumentation(*app_->Jnlst(),categories);
  }

  int numcols = getNumCols();
  if(obj_)
    delete [] obj_;
  obj_ = new double[numcols];
  CoinFillN(obj_,numcols,1.);
  setStrParam(OsiProbName, std::string(filename[1]));
  extractInterfaceParams();
  hasBeenOptimized_ = false;
  feasibilityProblem_ = new Ipopt::TNLP2FPNLP
      (Ipopt::SmartPtr<Ipopt::TNLP>(Ipopt::GetRawPtr(problem_)));
}

/** write ampl solution file */
void
BonminAmplInterface::writeAmplSolFile(std::string message,const double * primalSol,const double * dualSol)
{
  TMINLP * tminlp = GetRawPtr(tminlp_);
  AmplTMINLP * ampl_tminlp = dynamic_cast<AmplTMINLP *> (tminlp);
  if(ampl_tminlp)
    ampl_tminlp->write_solution(message,primalSol,dualSol);
  else
    std::cerr<<"Errot can not write .sol file for non ampl problem"<<std::endl;
}



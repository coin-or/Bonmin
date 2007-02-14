// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006 
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University
// Date:
// 06/29/2006

// Driver for Feasibility pump


#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <iomanip>
#include <fstream>

#include "CoinTime.hpp"

#include "BonAmplInterface.hpp"
#include "BonCbcParam.hpp"
#include "BonCbc.hpp"
#include "BonAmplTMINLP.hpp"
#include "AmplTNLP.hpp"
#include "FP.hpp"
#include "BonIpoptSolver.hpp"

void register_ALL_options( SmartPtr<RegisteredOptions> roptions );
void set_ipopt_minlp_default(SmartPtr<OptionsList> Option);



class AmplFP : public Bonmin::AmplTMINLP
{
  public:
  virtual void fillApplicationOptions2(Ipopt::AmplOptionsList* amplOptList)
{
    
    amplOptList->AddAmplOption("FP.Algo","FP.Algo",
                               AmplOptionsList::String_Option,
                               "specify minimal investment in share");
    amplOptList->AddAmplOption("FP.time_limit","FP.time_limit",
                               AmplOptionsList::Number_Option,
                               "Give time limit");
}

AmplFP(const SmartPtr<const Journalist>& jnlst, 
           const SmartPtr<OptionsList> options,
           char**& argv, 
           AmplSuffixHandler* suffix_handler /*=NULL*/,
           const std::string& appName,
           std::string* nl_file_content /* = NULL */):AmplTMINLP()
           {
               Initialize2(jnlst, options, argv, suffix_handler,
               appName, nl_file_content);
          }
virtual void
Initialize2(const SmartPtr<const Journalist>& jnlst, 
           const SmartPtr<OptionsList> options,
           char**& argv, 
           AmplSuffixHandler* suffix_handler =NULL,
           const std::string& appName = "fp",
           std::string* nl_file_content  = NULL )
           {
  SmartPtr<AmplOptionsList> ampl_options_list = new AmplOptionsList();
  fillAmplOptionList(GetRawPtr(ampl_options_list));
  fillApplicationOptions2(GetRawPtr(ampl_options_list) );
  std::string options_id = appName + "_options";
  ampl_tnlp_ = new AmplTNLP(jnlst, options, argv, suffix_handler, true,
                            ampl_options_list, options_id.c_str(),
                            appName.c_str(), appName.c_str(), nl_file_content);
           }
           public:
           AmplFP():AmplTMINLP(){}
           virtual AmplTMINLP * createEmpty(){return new AmplFP;}
};

class FPInterface : public Bonmin::AmplInterface
{
public:
  FPInterface(char **& argv):Bonmin::AmplInterface()
  {
    readAmplNlFile2(argv);
  }
  FPInterface(const FPInterface &other):
  Bonmin::AmplInterface(other)
  {}
  virtual OsiSolverInterface * clone()
  { return new FPInterface(*this);} 
protected:
  virtual std::string appName() {return "FP";}
  virtual void registerApplicationOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("FP options");
    roptions->AddStringOption4("Algo","type of FP algorithm",
                               "FP",
                               "FP", "Stand-alone FP",
                               "iFP", "iterated FP",
                               "OA", "Classical OA",
                               "eOA","Enhanced OA"
                               "");
  }
  void 
  readAmplNlFile2(char**& filename
                  )
  {
    app_ = new Bonmin::IpoptSolver();
    SmartPtr<RegisteredOptions> roptions = app_->RegOptions();
    register_ALL_options(roptions);
    registerApplicationOptions(roptions);
    
    // Call initalize to open output
    app_->Initialize("");
    // Read the bonmin.opt input file
    app_->Initialize("FP.opt");
    
    char * pbName = new char[strlen(filename[1])+1];
    strcpy(pbName, filename[1]);
    
    setStrParam(OsiProbName,std::string(pbName));
    delete [] pbName;
    
    
    if(!IsValid(tminlp_)) {
      amplTminlp_ = new AmplFP(ConstPtr(app_->Jnlst()), app_->Options(), filename,
                               NULL, appName() , NULL);
      tminlp_ = GetRawPtr(amplTminlp_);
    }
    else {
      Bonmin::AmplTMINLP * amplTMINLP = dynamic_cast<Bonmin::AmplTMINLP *> (GetRawPtr(tminlp_));
      if(amplTMINLP) {
        Bonmin::AmplTMINLP * newAmpl = amplTMINLP->createEmpty();
        newAmpl->Initialize(ConstPtr(app_->Jnlst()), app_->Options(), filename,
                            NULL, appName() , NULL);
        amplTminlp_ = newAmpl;
        tminlp_ = GetRawPtr(amplTminlp_);
      }
      else {
        amplTminlp_ = new AmplFP(ConstPtr(app_->Jnlst()), app_->Options(), filename,
                                     NULL, appName() , NULL);
        tminlp_ = GetRawPtr(amplTminlp_);
      }
    }
    problem_ = new Bonmin::TMINLP2TNLP(tminlp_);
    
    bool print_options_documentation;
    app_->Options()->GetBoolValue("print_options_documentation",
                                  print_options_documentation, "");
    if (print_options_documentation) {
      std::list<std::string> categories;
      categories.push_back("FP options");
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
    feasibilityProblem_ = new Bonmin::TNLP2FPNLP
      (Ipopt::SmartPtr<Ipopt::TNLP>(Ipopt::GetRawPtr(problem_)));
  }
};

int iteratedFP(Bonmin::AmplInterface &nlpSolver, 
               bool standAlone, 
               double *&solution);

int enhancedOA(Bonmin::AmplInterface &nlpSolver, 
               bool doFP, 
               double *&solution);

int main (int argc, char *argv[])
{
  using namespace Ipopt;
  
  FPInterface nlpSolver(argv);
  
                           
  //Set up done, now let's branch and bound
  //double time1 = CoinCpuTime();
  try {
    Ipopt::SmartPtr<Ipopt::OptionsList> Options = nlpSolver.retrieve_options();
    
    int algo;
    
    Options->GetEnumValue("Algo", algo, "FP.");
    Options->GetNumericValue("time_limit",params.maxTime_,"FP.");
    double * solution = NULL;
    if(algo==0)
      iteratedFP(nlpSolver,1, solution);
    else if (algo==1)
      iteratedFP(nlpSolver,0, solution);
    else if(algo==2)
      enhancedOA(nlpSolver,0, solution);
    else if(algo==3)
      enhancedOA(nlpSolver,1, solution);
    
    std::string message;
    if(solution==NULL) message="No solution";
      else message="Solution found";
    nlpSolver.writeAmplSolFile(message,solution,NULL);
    
  }
  catch(Bonmin::TNLPSolver::UnsolvedError *E) {
     E->printError(std::cerr);
    //There has been a failure to solve a problem with Ipopt.
    //And we will output file with information on what has been changed in the problem to make it fail.
    //Now depending on what algorithm has been called (B-BB or other) the failed problem may be at different place.
    //    const OsiSolverInterface &si1 = (algo > 0) ? nlpSolver : *model.solver();

  }
  catch(Bonmin::OsiTMINLPInterface::SimpleError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
    <<std::endl
    <<E.message()<<std::endl;
  }
  catch(CoinError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
    <<std::endl
    <<E.message()<<std::endl;
  }
  catch(...) {
    std::string pbName;
    
    nlpSolver.getStrParam(OsiProbName, pbName);

    std::cerr<<pbName<<" unrecognized excpetion"<<std::endl;
    std::cerr<<pbName<<"\t Finished \t exception"<<std::endl;
    throw;
  }
  
  return 0;
}



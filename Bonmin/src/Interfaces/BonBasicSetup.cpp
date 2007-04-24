// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 04/12/2007

#include "BonBasicSetup.hpp"
namespace Bonmin{
  
  BasicSetup::BasicSetup():
  journalist_(NULL),
  options_(NULL),
  roptions_(NULL)
{
  options_ = new Ipopt::OptionsList();
  
  journalist_= new Ipopt::Journalist();
  roptions_ = new Ipopt::RegisteredOptions();
  
  try{
    Ipopt::SmartPtr<Ipopt::Journal> stdout_journal =
    journalist_->AddFileJournal("console", "stdout", Ipopt::J_ITERSUMMARY);
    
    options_->SetJournalist(journalist_);
    options_->SetRegisteredOptions(roptions_);
  }
  catch (Ipopt::IpoptException &E){
    E.ReportException(*journalist_);
    throw E;
  }
  catch(std::bad_alloc){
    journalist_->Printf(Ipopt::J_ERROR, Ipopt::J_MAIN, "\n Not enough memory .... EXIT\n");
    throw -1;
  }
  catch(...){
    Ipopt::IpoptException E("Uncaught exception in FilterSolver::FilterSolver()",
                            "BonFilterSolver.cpp",-1);
    throw E;
  }
}

void BasicSetup::registerOptions(Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions){
  roptions->SetRegisteringCategory("Bonmin algotihm and solver choice");
  roptions->AddStringOption5("algorithm",
                             "Choice of the algorithm.",
                             "B-Hyb",
                             "B-BB","simple branch-and-bound algorithm,",
                             "B-OA","OA Decomposition algorithm,",
                             "B-QG","Quesada and Grossmann branch-and-cut algorithm,",
                             "B-Hyb","hybrid outer approximation based branch-and-cut.",
                             "B-Couenne","Branch-and-bound using Couenne convexifier (not available through Bonmin excutable)",
                             "This will preset default values for most options of bonmin but depending on which algorithm "
                             "some of these can be changed.");
  
  roptions->AddStringOption2("nlp_solver",
                             "Choice of the solver for local optima of continuous nlp's",
                             "Ipopt",
                             "Ipopt", "Interior Point OPTimizer (https://projects.coin-or.org/Ipopt)",
                             "filterSQP", "Sequential quadratic programming trust region algorithm (http://www-unix.mcs.anl.gov/~leyffer/solvers.html)",
                             "");

}

void 
BasicSetup::Initialize(const std::string& optFile)
{
  std::ifstream is;
  if (optFile != "") {
    try {
      is.open(optFile.c_str());
    }
    catch(std::bad_alloc) {
      journalist_->Printf(Ipopt::J_SUMMARY, Ipopt::J_MAIN, "\nEXIT: Not enough memory.\n");
      throw -1;
    }
    catch(...) {
      Ipopt::IpoptException E("Unknown Exception caught in ipopt", "Unknown File", -1);
      E.ReportException(*journalist_);
      throw -1;
    }
  }
  Initialize(is);
  if (is) {
    is.close();
  }
}

void 
BasicSetup::Initialize(std::istream& is)
{
  if(is.good()){
    options_->ReadFromStream(*journalist_, is);
  }
  bool print_options_documentation;
  options_->GetBoolValue("print_options_documentation",
                                print_options_documentation, "");
  if (print_options_documentation) {
    std::list<std::string> categories;
    categories.push_back("bonmin branch-and-bound options");
    categories.push_back("bonmin options for robustness");
    categories.push_back("bonmin options for non-convex problems");
    categories.push_back("bonmin options : B-Hyb specific options");
    //    roptions->OutputLatexOptionDocumentation2(*app_->Jnlst(),categories);
    roptions_->OutputOptionDocumentation(*(journalist_),categories);
  }
  
}
BaseOptions::Solver BasicSetup::getSolver(){
  int ival;
  options_->GetEnumValue("nlp_solver", ival,"bonmin.");
  return BaseOptions::Solver(ival);  
}

BaseOptions::Algorithm BasicSetup::getAlgorithm(){
  int ival;
  options_->GetEnumValue("algorithm", ival,"bonmin.");
  return BaseOptions::Algorithm(ival);  
}

void 
BasicSetup::InitializeFromLongString(std::string & string){
  std::stringstream is(string.c_str());
  Initialize(is);
}

/** May print the documentation if user asked.*/
void 
BasicSetup::mayPrintDoc(){
  bool print_options_documentation;
  options_->GetBoolValue("print_options_documentation",
                               print_options_documentation, "");
  if (print_options_documentation) {
    std::list<std::string> categories;
    categories.push_back("bonmin branch-and-bound options");
    categories.push_back("bonmin options for robustness");
    categories.push_back("bonmin options for non-convex problems");
    categories.push_back("bonmin options : B-Hyb specific options");
    //    roptions->OutputLatexOptionDocumentation2(*app_->Jnlst(),categories);
    roptions_->OutputOptionDocumentation(*journalist_,categories);
  }
}

}/** End namespace Bonmin. */


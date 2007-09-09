// (C) Copyright International Business Machines Corporation 2007
//
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines,
//
// Date : 27/08/2007

#include "BonRegisteredOptions.hpp"
#include "IpSmartPtr.hpp"
#include <sstream>

namespace Bonmin{
  struct optionsCmp{
  bool operator()( Ipopt::RegisteredOption * a,
                   Ipopt::RegisteredOption * b){
    if(a->RegisteringCategory() == b->RegisteringCategory()){
       return a->Name() < b->Name();
    }
    return a->RegisteringCategory() < b->RegisteringCategory(); 
  }
  };


  static std::string makeLatex(const std::string &s){
    std::string ret_val;
    std::string::const_iterator i = s.begin();
    for(; i != s.end() ; i++){
       switch (*i) {
         case '_':
         case '-':
           ret_val +='\\';
         default:
           ret_val += *i;
       }
    }
   return ret_val;
  }

  static std::string makeLatex(double value){
    std::string ret_val = "$";
    std::stringstream s_val;
    s_val<<value;
    
    unsigned int i = s_val.str().find_first_of('e');
    if(i != s_val.str().size()){
      ret_val += s_val.str().substr(0,i-1);
      ret_val += " \\cdot 10^{";
      ret_val += s_val.str().substr(i+1);
      ret_val += '}';
    }
    else ret_val += s_val.str();
    ret_val += '$';
    return ret_val;
  }
  static std::string makeLatex(int value){
    std::string ret_val = "$";
    std::stringstream s_val;
    s_val<<value;
    ret_val += s_val.str();
    ret_val += "$";
    return ret_val;
  }

  static char OptionType2Char(const Ipopt::RegisteredOptionType &T){
    switch(T){
      case Ipopt::OT_Number: return 'F';
      case Ipopt::OT_Integer: return 'I';
      case Ipopt::OT_String: return 'S';
      case Ipopt::OT_Unknown: 
      default: return 'U';
    }
  }

  static std::string defaultAsString(Ipopt::SmartPtr< Ipopt::RegisteredOption > opt){
    Ipopt::RegisteredOptionType T = opt->Type();
    switch(T){
      case Ipopt::OT_Number: return makeLatex(opt->DefaultNumber());
      case Ipopt::OT_Integer: return makeLatex(opt->DefaultInteger());
      case Ipopt::OT_String: return makeLatex(opt->DefaultString());
      case Ipopt::OT_Unknown: 
      default:
         return "Unknown type of option";
    }
  }
  /** Output Latex table of options.*/
  void 
  RegisteredOptions::writeLatexOptionsTable(std::ostream &of, ExtraCategoriesInfo which){
  std::map<std::string, Ipopt::SmartPtr<Ipopt::RegisteredOption> > 
           registered_options = RegisteredOptionsList();

  //Print table header
  of<<"\\begin{threeparttable}"<<std::endl
    <<"\\caption{\\label{tab:options} "<<std::endl
    <<"List of options and compatibility with the different algorithms."<<std::endl
    <<"}"<<std::endl
    <<"\\begin{tabular}{|l|r|r|r|r|r|r|}"
    <<"\\hline"<<std::endl
    <<"Option & type &  default & {\\tt B-BB} & {\\tt B-OA} & {\\tt B-QG} & {\\tt B-Hyb} \\\\"<<std::endl
    <<"\\hline"<<std::endl;

  //sort options by categories and alphabetical order
  std::list< Ipopt::RegisteredOption * > sortedOptions;
 
  for(std::map<std::string, Ipopt::SmartPtr<Ipopt::RegisteredOption > >::iterator i = 
          registered_options.begin(); i != registered_options.end() ; i++){
     if(categoriesInfo(i->second->RegisteringCategory()) == which)
     sortedOptions.push_back(GetRawPtr(i->second));
     }

   sortedOptions.sort(optionsCmp());
   std::string registeringCategory = "";
   for(std::list< Ipopt::RegisteredOption * >::iterator i = sortedOptions.begin();
       i != sortedOptions.end() ; i++)
   {
     if((*i)->RegisteringCategory() != registeringCategory){
     registeringCategory = (*i)->RegisteringCategory();
     of<<"\\hline"<<std::endl
       <<"\\multicolumn{1}{|c}{} & \\multicolumn{6}{l|}{"
       <<registeringCategory<<"}\\\\"<<std::endl
       <<"\\hline"<<std::endl;
     }
     
     of<<makeLatex((*i)->Name())<<"& "<<OptionType2Char((*i)->Type())<<"& "
       <<defaultAsString(*i)
       <<"& "<<( (isValidForBBB((*i)->Name()))? '+' : '-' )
       <<"& "<<( (isValidForBOA((*i)->Name()))? '+' : '-' )
       <<"& "<<( (isValidForBQG((*i)->Name()))? '+' : '-' )
       <<"& "<<( (isValidForHybrid((*i)->Name()))? '+' : '-' )
       <<"\\\\"<<std::endl;
   }
   //Print table end
  of<<"\\hline"<<std::endl
    <<"\\end{tabular}"<<std::endl
    <<"\\begin{tablenotes}"<<std::endl
    <<"\\item $\\strut^*$ option is available"<<std::endl
    <<"        for MILP subsolver (it is only passed if the {\\tt milp\\_subsolver} optio"<<std::endl
    <<"        see Subsection \\ref{sec:milp_opt})."<<std::endl
    <<"       \\item $\\strut^1$ disabled for stability reasons."<<std::endl
    <<"\\end{tablenotes}"<<std::endl
    <<"\\end{threeparttable} "<<std::endl;
  }

  /** choose options.*/
  void 
  RegisteredOptions::chooseOptions(ExtraCategoriesInfo which,
                                            std::list< Ipopt::RegisteredOption * >& sortedOptions)
  {
  std::map<std::string, Ipopt::SmartPtr<Ipopt::RegisteredOption> > 
           registered_options = RegisteredOptionsList();

  for(std::map<std::string, Ipopt::SmartPtr<Ipopt::RegisteredOption > >::iterator i = 
          registered_options.begin(); i != registered_options.end() ; i++){
     if(categoriesInfo(i->second->RegisteringCategory()) == which)
     sortedOptions.push_back(GetRawPtr(i->second));
     }
   sortedOptions.sort(optionsCmp());
  }
  /** Output html table of options.*/
  void 
  RegisteredOptions::writeHtmlOptionsTable(std::ostream &of, ExtraCategoriesInfo which){

  //Print table header
  of<<"<table border=\"1\">"<<std::endl;
  //sort options by categories and alphabetical order
  std::list< Ipopt::RegisteredOption * > sortedOptions;
  chooseOptions(which, sortedOptions); 
  writeHtmlOptionsTable(of, sortedOptions);
  }

   
   /** Output html table of options.*/
   void 
   RegisteredOptions::writeHtmlOptionsTable(std::ostream &os, std::list<Ipopt::RegisteredOption *> &options)
   {
   std::string registeringCategory = "";
   for(std::list< Ipopt::RegisteredOption * >::iterator i = options.begin();
       i != options.end() ; i++)
   {
     if((*i)->RegisteringCategory() != registeringCategory){
     registeringCategory = (*i)->RegisteringCategory();
     os<<"<tr>"
       <<"   <th colspan=7>"
       <<registeringCategory<<"</th>"<<std::endl
       <<"</tr>"<<std::endl;
     }
     
     os<<"<tr>"<<std::endl
       <<"<td> <a href=\"#"<<((*i)->Name())<<"\">"
       <<((*i)->Name())<<"</a> </td>"<<std::endl
       <<"<td>"<<OptionType2Char((*i)->Type())<<"</td>"<<std::endl
       <<"<td>"<<defaultAsString(*i)<<"</td>"<<std::endl
       <<"<td> "<<( (isValidForBBB((*i)->Name()))? '+' : '-' )<<"</td>"<<std::endl
       <<"<td>"<<( (isValidForBOA((*i)->Name()))? '+' : '-' )<<"</td>"<<std::endl
       <<"<td>"<<( (isValidForBQG((*i)->Name()))? '+' : '-' )<<"</td>"<<std::endl
       <<"<td>"<<( (isValidForHybrid((*i)->Name()))? '+' : '-' )<<"</td>"<<std::endl
       <<"</tr>"<<std::endl;
   }
   //Print table end
  os<<"</tr>"<<std::endl
    <<"</table>"<<std::endl;
  }

//Copiedand modified from Ipopt/src/Common adapt to output to an ostream
  static void 
  OutputLatexDescription(std::ostream & os, const Ipopt::RegisteredOption opt)
  {
    std::string latex_name = makeLatex(opt.Name());
    std::string latex_desc = makeLatex(opt.ShortDescription());
    os<<"\\paragraph{"<<latex_name<<"}\n"
      <<"\\label{"<<latex_name<<"}\n"
      <<latex_desc<<"\\\\"<<std::endl;

    os<<makeLatex(opt.LongDescription());

    Ipopt::RegisteredOptionType type = opt.Type();
    if (type == Ipopt::OT_Number) {
      os<<" The valid range for this real option is \n$";
      if (opt.HasLower()) {
         os<<makeLatex(opt.LowerNumber());
      }
      else {
         os<<"-\\infty$";
      }

      if (opt.HasLower() && !opt.LowerStrict()) {
         os<<" \\le ";
      }
      else {
        os<<" <  ";
      }

      os<<"{\\tt "<<latex_name<<"}";

      if (opt.HasUpper() && !opt.UpperStrict()) {
        os<<" \\le ";
      }
      else {
        os<<" <  ";
      }

      if (opt.HasUpper()) {
        os<<makeLatex(opt.UpperNumber());
      }
      else {
        os<<"\\infty";
      }

      os<<"$\nand its default value is $"<<
          opt.DefaultNumber()<<"$.\n\n";
    }
    else if (type == Ipopt::OT_Integer) {
      os<<" The valid range for this integer option is\n$";
      if (opt.HasLower()) {
        os<<opt.LowerInteger()<<" \\le ";
      }
      else {
        os<<"{\\tt - \\infty}";
      }

      os<< "{\\tt "<<latex_name<<"}";

      if (opt.HasUpper()) {
        os<<" \\le "<<opt.UpperInteger();
      }
      else {
        os<<" <  {\\tt \\infty}";
      }

      os<<"$\nand its default value is $"
        <<opt.DefaultInteger()<<"$.\n\n";
    }
    else if (type == Ipopt::OT_String) {
      std::string buff;
      os<<"\nThe default value for this string option is {\\tt "
         <<opt.DefaultString()<<"}.\n";

      os<<"\\\\ \nPossible values:\n";
      os<<"\\begin{itemize}\n";
#if TOFIX
      for (std::vector<string_entry>::const_iterator
           i = valid_strings_.begin();
           i != valid_strings_.end(); i++) {
        std::string latex_value;
        MakeValidLatexString((*i).value_, latex_value);
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "   \\item %s: ",
                     latex_value.c_str());

        std::string latex_desc;
        MakeValidLatexString((*i).description_, latex_desc);
        jnlst.PrintStringOverLines(J_SUMMARY, J_DOCUMENTATION, 0, 48,
                                   latex_desc.c_str());
        jnlst.Printf(J_SUMMARY, J_DOCUMENTATION, "\n");
      }
#endif
      std::cout<<"\\end{itemize}\n"<<std::endl;
    }
   std::cout<<std::endl;
  }
   /** Output Latex/Html ooptions documentation.*/
   void 
   RegisteredOptions::writeLatexHtmlDoc(std::ostream &os, ExtraCategoriesInfo which){
#ifdef FIXME
      std::list< Ipopt::RegisteredOptions * > options;
      ExtraCategoriesInfo which = Bonmin
      chooseOptions(which, options);
      os<<"\\htmlonly{"
      os<<"\\begin{rawhtml}"<<std::endl;
      writeHtmlOptionsTable(os, options)
      os<<"\\end{rawhtml} }"<<std::endl;

      std::string registeringCategory = "";
      for(std::list< Ipopt::RegisteredOption * >::iterator i = options.begin();
           i != options.end() ; i++)
       {
          if((*i)->RegisteringCategory() != registeringCategory){
           registeringCategory = (*i)->RegisteringCategory();
             os<<"\\subsection{"<<registeringCategory<<"}"<<std::endl;      
             os<<"\\label{sec:"<<registeringCategory<<"}"<<std::endl;
           }
       
       //MISSING CODE     
       }
#endif
    }

}/*Ends bonmin namespace.*/

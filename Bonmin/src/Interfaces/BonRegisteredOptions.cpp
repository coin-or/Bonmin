// (C) Copyright International Business Machines Corporation 2007
//
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, International Business Machines,
//
// Date : 27/08/2007

#include "BonRegisteredOptions.hpp"
#include <IpSmartPtr.hpp>
#include <sstream>
#include <climits>
#include <cfloat>

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

  static std::string makeSpaceLess(const std::string &s){
    std::string ret_val;
    std::string::const_iterator i = s.begin();
    for(; i != s.end() ; i++){
       switch (*i) {
         case ' ':
         case '\t':
         case '_':
           //ret_val +='_';
           break;
         default:
           ret_val += *i;
       }
    }
   return ret_val;
  }

#if 0
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
#endif

  static std::string makeString(int value){
    std::string ret_val;
    if(value >= INT_MAX){
     ret_val="INT_MAX";}
    else if(value <= - INT_MAX){
      ret_val="-INT_MAX";}
    else{
      std::stringstream s_val;
      s_val<<value;
      ret_val = s_val.str();
    }
    return ret_val;
  }

  static std::string makeString(double value){
    std::string ret_val;
    if(value >= DBL_MAX){
     ret_val="DBL_MAX";}
    else if(value <= - DBL_MAX){
      ret_val="-DBL_MAX";}
    else{
      std::stringstream s_val;
      s_val<<value;
      ret_val = s_val.str();
    }
    return ret_val;
  }

 static std::string makeNumber(std::string value){
   if(value == "DBL_MAX") {
     std::stringstream s_val;
     s_val<<DBL_MAX;
     return s_val.str();
   }
   if(value == "-DBL_MAX") {
     std::stringstream s_val;
     s_val<<-DBL_MAX;
     return s_val.str();
   }
   if(value == "INT_MAX") {
     std::stringstream s_val;
     s_val<<INT_MAX;
     return s_val.str();
   }
   if(value == "-INT_MAX") {
     std::stringstream s_val;
     s_val<<-INT_MAX;
     return s_val.str();
   }
   return value;
 }
#if 0
  static std::string makeLatex(int value){
    std::string ret_val = "$";
    std::stringstream s_val;
    s_val<<value;
    ret_val += s_val.str();
    ret_val += "$";
    return ret_val;
  }
#endif

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
      case Ipopt::OT_Number: return makeString(opt->DefaultNumber());
      case Ipopt::OT_Integer: return makeString(opt->DefaultInteger());
      case Ipopt::OT_String: return (opt->DefaultString());
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
  //of<<"\\begin{threeparttable}"<<std::endl
  of<<"\\topcaption{\\label{tab:options} "<<std::endl
    <<"List of options and compatibility with the different algorithms."<<std::endl
    <<"}"<<std::endl;
  of<<"\\tablehead{\\hline "<<std::endl
    <<"Option & type & ";
  //of<<" default & ";
  of<<"{\\tt B-BB} & {\\tt B-OA} & {\\tt B-QG} & {\\tt B-Hyb} & {\\tt B-Ecp} & {\\tt B-iFP} & {\\tt Cbc\\_Par} \\\\"<<std::endl
    <<"\\hline"<<std::endl
    <<"\\hline}"<<std::endl;
  of<<"\\tabletail{\\hline \\multicolumn{9}{|c|}{continued on next page}\\\\"
    <<"\\hline}"<<std::endl; 
  of<<"\\tablelasttail{\\hline}"<<std::endl;
  of<<"{\\footnotesize"<<std::endl;
  of<<"\\begin{xtabular}{@{}|@{\\;}l@{\\;}|@{\\;}r@{\\;}|@{\\;}c@{\\;}|@{\\;}c@{\\;}|@{\\;}c@{\\;}|@{\\;}c@{\\;}|@{\\;}c@{\\;}|@{\\;}c@{\\;}|@{\\;}c@{\\;}|@{}}"<<std::endl;

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
       <<"\\multicolumn{1}{|c}{} & \\multicolumn{8}{l|}{"
       <<registeringCategory<<"}\\\\"<<std::endl
       <<"\\hline"<<std::endl;
     }
     
     of<<makeLatex((*i)->Name())<<"& "<<OptionType2Char((*i)->Type())//<<"& "
       //<<makeLatex(defaultAsString(*i))
       <<"& "<<( (isValidForBBB((*i)->Name()))? "$\\surd$" : "-" )
       <<"& "<<( (isValidForBOA((*i)->Name()))? "$\\surd$" : "-" )
       <<"& "<<( (isValidForBQG((*i)->Name()))? "$\\surd$" : "-" )
       <<"& "<<( (isValidForHybrid((*i)->Name()))? "$\\surd$" : "-" )
       <<"& "<<( (isValidForBEcp((*i)->Name()))? "$\\surd$" : "-" )
       <<"& "<<( (isValidForBiFP((*i)->Name()))? "$\\surd$" : "-" )
       <<"& "<<( (isValidForCbc((*i)->Name()))? "$\\surd$" : "-" )
       <<"\\\\"<<std::endl;
   }
   //Print table end
  of<<"\\hline"<<std::endl
    <<"\\end{xtabular}"<<std::endl;
  of<<"}"<<std::endl;
#if 0
  of<<"\\begin{tablenotes}"<<std::endl
    <<"\\item $\\strut^*$ option is available"<<std::endl
    <<"        for MILP subsolver (it is only passed if the {\\tt milp\\_subsolver} optio"<<std::endl
    <<"        see Subsection \\ref{sec:milp_opt})."<<std::endl
    <<"       \\item $\\strut^1$ disabled for stability reasons."<<std::endl
    <<"\\end{tablenotes}"<<std::endl
    <<"\\end{threeparttable} "<<std::endl;
#endif
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
   os<<"<table border=\"1\">"<<std::endl;
   os<<"<tr>"<<std::endl;
   os<<"<td>Option </td>"<<std::endl;
   os<<"<td> type </td>"<<std::endl;
   //os<<"<td> default </td>"<<std::endl;
   os<<"<td> B-BB</td>"<<std::endl;
   os<<"<td> B-OA</td>"<<std::endl;
   os<<"<td> B-QG</td>"<<std::endl;
   os<<"<td> B-Hyb</td>"<<std::endl;
   os<<"</tr>"<<std::endl;
   std::string registeringCategory = "";
   for(std::list< Ipopt::RegisteredOption * >::iterator i = options.begin();
       i != options.end() ; i++)
   {
     if((*i)->RegisteringCategory() != registeringCategory){
     registeringCategory = (*i)->RegisteringCategory();
     os<<"<tr>"
       <<"   <th colspan=9>"
       <<" <a href=\"#sec:"<<makeSpaceLess(registeringCategory)<<"\">"
       <<registeringCategory<<"</a> </th>"<<std::endl
       <<"</tr>"<<std::endl;
     }
     
     os<<"<tr>"<<std::endl
       <<"<td>"<<((*i)->Name())<<"</td>"<<std::endl
       <<"<td>"<<OptionType2Char((*i)->Type())<<"</td>"<<std::endl
       //<<"<td>"<<defaultAsString(*i)<<"</td>"<<std::endl
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

   /** Output Latex/Html options documentation.*/
   void 
   RegisteredOptions::writeLatexHtmlDoc(std::ostream &os, ExtraCategoriesInfo which){
      std::list< Ipopt::RegisteredOption * > options;
      chooseOptions(which, options);
      os<<"\\latexhtml{}{"<<std::endl;
      os<<"\\HCode{"<<std::endl;
      writeHtmlOptionsTable(os, options);
      os<<"}\n}"<<std::endl;

      //Create journalist to write to os
      Ipopt::Journalist jnlst;
      Ipopt::SmartPtr<Ipopt::StreamJournal> J = new Ipopt::StreamJournal("options_journal", Ipopt::J_ALL);
      J->SetOutputStream(&os);
      J->SetPrintLevel(Ipopt::J_DOCUMENTATION, Ipopt::J_SUMMARY);
      jnlst.AddJournal(GetRawPtr(J));

      std::string registeringCategory = "";
      for(std::list< Ipopt::RegisteredOption * >::iterator i = options.begin();
           i != options.end() ; i++)
       {
          if((*i)->RegisteringCategory() != registeringCategory){
           registeringCategory = (*i)->RegisteringCategory();
             os<<"\\subsection{"<<registeringCategory<<"}"<<std::endl;      
             os<<"\\label{sec:"<<makeSpaceLess(registeringCategory)<<"}"<<std::endl;
             os<<"\\htmlanchor{sec:"<<makeSpaceLess(registeringCategory)<<"}"<<std::endl;
           }
       
           (*i)->OutputLatexDescription(jnlst);
       }
    }

  /** Ouptut a bonmin.opt file with options default values and short descritpions.*/
  void
  RegisteredOptions::writeBonminOpt(std::ostream &os, ExtraCategoriesInfo which){
      std::list< Ipopt::RegisteredOption * > options;
      chooseOptions(which, options);

      //Create journalist to write to os
      Ipopt::Journalist jnlst;
      Ipopt::SmartPtr<Ipopt::StreamJournal> J = new Ipopt::StreamJournal("options_journal", Ipopt::J_ALL);
      J->SetOutputStream(&os);
      J->SetPrintLevel(Ipopt::J_DOCUMENTATION, Ipopt::J_SUMMARY);
      jnlst.AddJournal(GetRawPtr(J));

      std::string registeringCategory = "";
      for(std::list< Ipopt::RegisteredOption * >::iterator i = options.begin();
           i != options.end() ; i++)
       {
          if((*i)->RegisteringCategory() != registeringCategory){
           registeringCategory = (*i)->RegisteringCategory();
             os<<std::endl<<"# registering category: "<<registeringCategory<<std::endl<<std::endl;
           }
           os<<"bonmin.";
           os.setf(std::ios::left);
           os.width(37);
           os<<(*i)->Name()<<" ";
           os.width(10);
           os<<makeNumber(defaultAsString(*i))<<"\t#";
           os<<(*i)->ShortDescription();
           os<<std::endl;            
       }
    }

}/*Ends bonmin namespace.*/

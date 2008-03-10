// (C) Copyright International Business Machines Corporation 2007
//
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 27/08/2007

#ifndef BonRegisteredOptions_H
#define BonRegisteredOptions_H

#include "IpRegOptions.hpp"
#include "IpException.hpp"
#include "IpTypes.hpp"
#include "CoinError.hpp"
#include <iostream>
namespace Bonmin {
/** Class to add a few more information to Ipopt::RegisteredOptions.
    In particular, it allows to store code to indicate in which algorithm
    option is available. It also allows to table summing up all the options
    both in LaTex and html.*/
class RegisteredOptions: public Ipopt::RegisteredOptions{
  public:
    enum ExtraOptInfosBits{
    validInHybrid=0/** Say that option is valid in Hybrid method.*/,
    validInQG/** Say that option is valid in Quesada Grossmann method.*/,
    validInOA/**Say that option is valid in outer approximation dec.*/,
    validInBBB/** Say that option is valid in the pure branch-and-bound.*/
   };
   enum ExtraCategoriesInfo{
    BonminCategory = 0/** Option category is for Bonmin.*/,
    IpoptCategory /** Option category for Ipopt.*/,
    FilterCategory /** Option category for FilterSqp.*/,
    BqpdCategory /** Option category for Bqpd.*/,
    CouenneCategory /** Option category for Couenne.*/
   };
    /** Standard constructor.*/
    RegisteredOptions():
       Ipopt::RegisteredOptions(){
    }

    /** Standard destructor.*/
    ~RegisteredOptions(){
    }

   //DECLARE_STD_EXCEPTION(OPTION_NOT_REGISTERED); 
   /** Set registering category with extra information.*/
   void SetRegisteringCategory (const std::string &registering_category,
                                const ExtraCategoriesInfo extra){
      Ipopt::RegisteredOptions::SetRegisteringCategory(registering_category);
      categoriesInfos_[registering_category] = extra;}
 
   /** throw if option does not exists.*/
   inline void optionExists(const std::string & option){
     if(!IsValid(GetOption(option))){
       std::string msg = "Try to access option: "+option;
       msg += "\n Option is not registered.\n";
       throw CoinError("Bonmin::RegisteredOption","optionExists",msg);
     }
   }

   /**Set extra information for option.*/
   inline void setOptionExtraInfo(const std::string & option, int code){
      optionExists(option);
      bonOptInfos_[option] = code;
   }

   /** Set that option is valid for hybrid.*/
   inline void optionValidForHybrid(const std::string &option){
      optionExists(option);
     bonOptInfos_[option] |= 1 << validInHybrid;}

   /** Set that option is valid for QuesadaGrossmann.*/
   inline void optionValidForBQG(const std::string &option){
     optionExists(option);
     bonOptInfos_[option] |= 1 << validInQG;}
   
   /** Set that option is valid for Outer approximation.*/
   inline void optionValidForBOA(const std::string &option){
     optionExists(option);
     bonOptInfos_[option] |= 1 << validInOA;}
   
   /** Set that option is valid for pure branch-and-bound.*/
   inline void optionValidForBBB(const std::string &option){
     optionExists(option);
     bonOptInfos_[option] |= 1 << validInBBB;}


    /** Say if option is valid for hybrid.*/
   inline bool isValidForHybrid(const std::string &option){
      optionExists(option);
     std::map<std::string, int>::iterator i = bonOptInfos_.find(option);
     if(i != bonOptInfos_.end()) 
     return (i->second) & (1 << validInHybrid);
     else return true;}

   /** Say if option is valid for QuesadaGrossmann.*/
   inline bool isValidForBQG(const std::string &option){
     optionExists(option);
     std::map<std::string, int>::iterator i = bonOptInfos_.find(option);
     if(i != bonOptInfos_.end()) 
       return (i->second) & (1 << validInQG);
     else return true;}
   
   /** Say if option is valid for Outer approximation.*/
   inline bool isValidForBOA(const std::string &option){
     optionExists(option);
     std::map<std::string, int>::iterator i = bonOptInfos_.find(option);
     if(i != bonOptInfos_.end()) 
     return (i->second) & (1 << validInOA);
     return true;}
   
   /** Say if option is valid for pure branch-and-bound.*/
   inline bool isValidForBBB(const std::string &option){
     optionExists(option);
     std::map<std::string, int>::iterator i = bonOptInfos_.find(option);
     if(i != bonOptInfos_.end()) 
     return (i->second) & (1 << validInBBB);
     return true;}
   /** Output Latex table of options.*/
   void writeLatexOptionsTable(std::ostream &of, ExtraCategoriesInfo which);

   /** Output html table of options.*/
   void writeHtmlOptionsTable(std::ostream &of, ExtraCategoriesInfo which);


   /** Output Latex/Html ooptions documentation.*/
   void writeLatexHtmlDoc(std::ostream &of, ExtraCategoriesInfo which);
   private:
   /** Output Latex table of options.*/
   void chooseOptions(ExtraCategoriesInfo which, std::list<Ipopt::RegisteredOption *> &options);
   /** Output html table of options.*/
   void writeHtmlOptionsTable(std::ostream & os, std::list<Ipopt::RegisteredOption *> &options);
   /** Store extra Informations on Bonmin options.*/
   std::map<std::string, int> bonOptInfos_;
   /** Store extra Informations on Registering categories
       (is bonmin, filterSqp...).*/
   std::map<std::string, ExtraCategoriesInfo> categoriesInfos_;
   /** Get categoriesInfo taking care of Ipopt.*/
   ExtraCategoriesInfo categoriesInfo(const std::string &s)
   {
      std::map<std::string, ExtraCategoriesInfo>::iterator i = categoriesInfos_.find(s);
      if(i == categoriesInfos_.end())
        return IpoptCategory;
      return i->second;
   }
};

}/*Ends namespace Bonmin.*/
#endif


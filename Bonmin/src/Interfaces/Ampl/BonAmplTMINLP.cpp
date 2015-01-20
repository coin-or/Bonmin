// (C) Copyright International Business Machines Corporation and
// Carnegie Mellon University 2004-2011
//
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Carl D. Laird, Carnegie Mellon University,
// Andreas Waechter, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 12/01/2004
#include "IpBlas.hpp"

#include <list>

#include "AmplTNLP.hpp"
#include "BonAmplTMINLP.hpp"
#include <iostream>
#include <fstream>

#include "CoinHelperFunctions.hpp"
#include "BonExitCodes.hpp"

using namespace Ipopt;

namespace Bonmin{
  void 
  RegisteredOptions::fillAmplOptionList(ExtraCategoriesInfo which, Ipopt::AmplOptionsList * amplOptList){
      std::list<int> test;
      std::list< Ipopt::RegisteredOption * > options;
      chooseOptions(which, options);
      for(std::list< Ipopt::RegisteredOption * >::iterator i = options.begin();
           i != options.end() ; i++)
       {
           std::string name = "bonmin.";
           name += (*i)->Name();
           Ipopt::RegisteredOptionType T = (*i)->Type();
           Ipopt::AmplOptionsList::AmplOptionType  type;
           switch(T){
             case Ipopt::OT_Number: type = Ipopt::AmplOptionsList::Number_Option;
                  break;
             case Ipopt::OT_Integer: type = Ipopt::AmplOptionsList::Integer_Option;
                  break;
             case Ipopt::OT_String: type = Ipopt::AmplOptionsList::String_Option;
                  break;
             case Ipopt::OT_Unknown:
             default:
                throw CoinError("RegisteredOptions","fillAmplOptionList","Unknown option type");
           }
           amplOptList->AddAmplOption(name, name, type, (*i)->ShortDescription());
       }
}
}
#include "asl.h"
#include "asl_pfgh.h"
#include "getstub.h"

namespace ampl_utils
{
  void sos_kludge(int nsos, int *sosbeg, double *sosref);
}
namespace Bonmin
{

  AmplTMINLP::AmplTMINLP()
      :
      TMINLP(),
      appName_(),
      upperBoundingObj_(-1),
      ampl_tnlp_(NULL),
      branch_(),
      sos_(),
      suffix_handler_(NULL),
      constraintsConvexities_(NULL),
      numberNonConvex_(0),
      nonConvexConstraintsAndRelaxations_(NULL),
      numberSimpleConcave_(0),
      simpleConcaves_(NULL),
      hasLinearObjective_(false)
  {}


  AmplTMINLP::AmplTMINLP(const SmartPtr<const Journalist>& jnlst,
      const SmartPtr<Bonmin::RegisteredOptions> roptions,
      const SmartPtr<OptionsList> options,
      char**& argv,
      AmplSuffixHandler* suffix_handler /*=NULL*/,
      const std::string& appName,
      std::string* nl_file_content /* = NULL */
                        )
      :
      TMINLP(),
      appName_(),
      upperBoundingObj_(-1),
      ampl_tnlp_(NULL),
      branch_(),
      sos_(),
      suffix_handler_(NULL),
      constraintsConvexities_(NULL),
      numberNonConvex_(0),
      nonConvexConstraintsAndRelaxations_(NULL),
      numberSimpleConcave_(0),
      simpleConcaves_(NULL),
      hasLinearObjective_(false)
  {
    Initialize(jnlst, roptions, options, argv, suffix_handler, appName, nl_file_content);
  }

  void
  AmplTMINLP::Initialize(const SmartPtr<const Journalist>& jnlst,
      const SmartPtr<Bonmin::RegisteredOptions> roptions,
      const SmartPtr<OptionsList> options,
      char**& argv,
      AmplSuffixHandler* suffix_handler /*=NULL*/,
      const std::string& appName,
      std::string* nl_file_content /* = NULL */
                        )
  {
    appName_ = appName;
    options->GetEnumValue("file_solution",writeAmplSolFile_,"bonmin.");
    jnlst_ = jnlst;

    if (suffix_handler==NULL)
      suffix_handler_ = suffix_handler = new AmplSuffixHandler();

    // Add the suffix handler for scaling
    suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Objective_Source, AmplSuffixHandler::Number_Type);

    // Modified for warm-start from AMPL
    suffix_handler->AddAvailableSuffix("ipopt_zL_out", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("ipopt_zU_out", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("ipopt_zL_in", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("ipopt_zU_in", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);

    // priority suffix
    suffix_handler->AddAvailableSuffix("priority", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
    suffix_handler->AddAvailableSuffix("direction", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("downPseudocost", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("upPseudocost", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);



    // sos suffixes
    suffix_handler->AddAvailableSuffix("ref", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("sos",AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
    suffix_handler->AddAvailableSuffix("sos",AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Index_Type);
    suffix_handler->AddAvailableSuffix("sosno",AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("sosref",AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("sstatus", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
    suffix_handler->AddAvailableSuffix("sstatus", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Index_Type);


    // For marking convex/nonconvex constraints
    suffix_handler->AddAvailableSuffix("non_conv",AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
    suffix_handler->AddAvailableSuffix("primary_var",AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Index_Type);

    // For marking onoff constraints and indicator variables
    suffix_handler->AddAvailableSuffix("onoff_c",AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Index_Type);
    suffix_handler->AddAvailableSuffix("onoff_v",AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);

    // For objectives
    suffix_handler->AddAvailableSuffix("UBObj", AmplSuffixHandler::Objective_Source, AmplSuffixHandler::Index_Type);


    // Perturbation radius
    suffix_handler->AddAvailableSuffix("perturb_radius",AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);

    SmartPtr<AmplOptionsList> ampl_options_list = new AmplOptionsList();
    roptions->fillAmplOptionList(RegisteredOptions::BonminCategory, GetRawPtr(ampl_options_list));
    roptions->fillAmplOptionList(RegisteredOptions::FilterCategory, GetRawPtr(ampl_options_list));
    roptions->fillAmplOptionList(RegisteredOptions::BqpdCategory, GetRawPtr(ampl_options_list));
    fillApplicationOptions(GetRawPtr(ampl_options_list) );
    std::string options_id = appName + "_options";
    ampl_tnlp_ = new AmplTNLP(jnlst, options, argv, suffix_handler, true,
        ampl_options_list, options_id.c_str(),
        appName.c_str(), appName.c_str(), nl_file_content);


    /* Read suffixes */
    read_obj_suffixes();
    read_priorities();
    read_convexities();
    read_onoff();
    read_sos();


    /* Determine if objective is linear.*/
    Index n_non_linear_b= 0;
    Index n_non_linear_bi= 0;
    Index n_non_linear_c= 0;
    Index n_non_linear_ci = 0;
    Index n_non_linear_o= 0;
    Index n_non_linear_oi = 0;
    Index n_binaries = 0;
    Index n_integers = 0;
    ampl_tnlp_->get_discrete_info(n_non_linear_b, n_non_linear_bi, n_non_linear_c,
        n_non_linear_ci, n_non_linear_o, n_non_linear_oi,
        n_binaries, n_integers);
    if (n_non_linear_b == 0 && n_non_linear_o == 0) {
      hasLinearObjective_ = true;
    }
  }

  AmplTMINLP::~AmplTMINLP()
  {
    delete [] constraintsConvexities_;
    delete [] simpleConcaves_;
    delete [] nonConvexConstraintsAndRelaxations_;
    delete ampl_tnlp_;
  }

  void
  AmplTMINLP::read_priorities()
  {
    int numcols, m, dummy1, dummy2;
    TNLP::IndexStyleEnum index_style;
    ampl_tnlp_->get_nlp_info(numcols, m, dummy1, dummy2, index_style);

    const AmplSuffixHandler * suffix_handler = GetRawPtr(suffix_handler_);

    const Index* pri = suffix_handler->GetIntegerSuffixValues("priority", AmplSuffixHandler::Variable_Source);
    const Index* brac = suffix_handler->GetIntegerSuffixValues("direction", AmplSuffixHandler::Variable_Source);
    const Number* upPs = suffix_handler->GetNumberSuffixValues("upPseudocost", AmplSuffixHandler::Variable_Source);
    const Number* dwPs = suffix_handler->GetNumberSuffixValues("downPseudocost", AmplSuffixHandler::Variable_Source);


    branch_.gutsOfDestructor();
    branch_.size = numcols;
    if (pri) {
      branch_.priorities = new int[numcols];
      for (int i = 0 ; i < numcols ; i++) {
        branch_.priorities [i] = -pri[i] + 9999;
      }
    }
    if (brac) {
      branch_.branchingDirections = CoinCopyOfArray(brac,numcols);
    }
    if (upPs && !dwPs) dwPs = upPs;
    else if (dwPs && !upPs) upPs = dwPs;

    if (upPs) {
      branch_.upPsCosts = CoinCopyOfArray(upPs,numcols);
    }
    if (dwPs) {
      branch_.downPsCosts = CoinCopyOfArray(dwPs,numcols);
    }

    const double* perturb_radius =
      suffix_handler->GetNumberSuffixValues("perturb_radius", AmplSuffixHandler::Variable_Source);
    perturb_info_.SetPerturbationArray(numcols, perturb_radius);
  }

  void
  AmplTMINLP::read_sos()
  {
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();
    int i = 0;//ASL_suf_sos_explict_free;
    int copri[2], **p_sospri;
    copri[0] = 0;
    copri[1] = 0;
    int * starts = NULL;
    int * indices = NULL;
    char * types = NULL;
    double * weights = NULL;
    int * priorities = NULL;
    p_sospri = &priorities;

    sos_.gutsOfDestructor();
    int m = n_con;
    sos_.num = suf_sos(i, &sos_.numNz, &types, p_sospri, copri,
        &starts, &indices, &weights);
    if(m != n_con){
      throw CoinError("number of constraints changed by suf_sos. Not supported.",
                       "read_sos","Bonmin::AmplTMINLP");
   }
    if (sos_.num) {
      //Copy sos information
      sos_.priorities = CoinCopyOfArray(priorities,sos_.num);
      sos_.starts = CoinCopyOfArray(starts, sos_.num + 1);
      sos_.indices = CoinCopyOfArray(indices, sos_.numNz);
      sos_.types = CoinCopyOfArray(types, sos_.num);
      sos_.weights = CoinCopyOfArray(weights, sos_.numNz);

      ampl_utils::sos_kludge(sos_.num, sos_.starts, sos_.weights);
      for (int ii=0;ii<sos_.num;ii++) {
        int ichar = sos_.types[ii] - '0';
        if (ichar != 1 && ichar != 2) {
          std::cerr<<"Unsuported type of sos constraint: "<<sos_.types[ii]<<std::endl;
          throw;
        }
        sos_.types[ii]= static_cast<char>(ichar);
      }
    }
  }

  void
  AmplTMINLP::read_obj_suffixes()
  {
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();
    DBG_ASSERT(asl);
    //Get the values
    if (n_obj < 2) return;

    const AmplSuffixHandler * suffix_handler = GetRawPtr(suffix_handler_);

    const Index* UBObj = suffix_handler->GetIntegerSuffixValues("UBObj", AmplSuffixHandler::Objective_Source);
    if (UBObj) {
      //get index of lower bounding objective
      int lowerBoundingObj = (UBObj[0] == 1)? 1:0;
      // Pass information to Ipopt
      ampl_tnlp_->set_active_objective(lowerBoundingObj);

      //get the index of upper bounding objective
      for (int i = 0; i < n_obj; i++) {
        if (UBObj[i]==1) {
          if (upperBoundingObj_ != -1) {
            jnlst_->Printf(J_ERROR, J_MAIN,
                "Too many objectives for upper-bounding");
          }
          upperBoundingObj_ = i;
        }
      }
    }
    else {
      ampl_tnlp_->set_active_objective(0);
    }
  }
  /** To store all data stored in the nonconvex suffixes.*/
  struct NonConvexSuff
  {
    /** Default constructor.*/
    NonConvexSuff():
        cIdx(-1),relIdx(-1),scXIdx(-1),scYIdx(-1)
    {}
    /** Copy constructor.*/
    NonConvexSuff(const NonConvexSuff& other):
        cIdx(other.cIdx), relIdx(other.relIdx),
        scXIdx(other.scXIdx), scYIdx(other.scYIdx)
    {}
    /** Index of the nonconvex constraint.*/
    int cIdx;
    /** Index of its relaxation.*/
    int relIdx;
    /** Index of variable x in a simple concave constraint of type y >= F(x).*/
    int scXIdx;
    /** Index of variable y in a simple concave constraint of type y >= F(x).*/
    int scYIdx;
  };
  void AmplTMINLP::read_convexities()
  {
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();
    DBG_ASSERT(asl);

    const AmplSuffixHandler * suffix_handler = GetRawPtr(suffix_handler_);
    const Index * id = suffix_handler->GetIntegerSuffixValues("non_conv", AmplSuffixHandler::Variable_Source);
    const Index * primary_var = suffix_handler->GetIntegerSuffixValues("primary_var", AmplSuffixHandler::Constraint_Source);


    if (primary_var!= NULL) {
      if (constraintsConvexities_ != NULL) {
        delete [] constraintsConvexities_;
      }
      constraintsConvexities_ = new TMINLP::Convexity[n_con];
      if (id == NULL) {
        std::cerr<<"Incorrect suffixes description in ampl model. n_conv's are not declared "<<std::endl;
        exit(ERROR_IN_AMPL_SUFFIXES);
      }
      int numberSimpleConcave = 0;
      std::map<int, int> id_map;

      for (int i = 0 ; i < n_var ; i++) {
        id_map[id[i]] = i;
      }


      for (int i = 0 ; i < n_con ; i++) {
        if (primary_var[i] != 0) {
          constraintsConvexities_[i] = TMINLP::SimpleConcave;
          numberSimpleConcave++;
        }
        else constraintsConvexities_[i] = TMINLP::Convex;
      }
      simpleConcaves_ = new SimpleConcaveConstraint[numberSimpleConcave];
      nonConvexConstraintsAndRelaxations_ = new MarkedNonConvex[numberSimpleConcave];
      numberSimpleConcave = 0;
      int * jCol = new int[n_var];
      for (int i = 0 ; i < n_con ; i++) {
        if (primary_var[i] != 0) {
          nonConvexConstraintsAndRelaxations_[numberSimpleConcave].cIdx = i;
          nonConvexConstraintsAndRelaxations_[numberSimpleConcave].cRelaxIdx = -1;
          simpleConcaves_[numberSimpleConcave].cIdx = i;
          simpleConcaves_[numberSimpleConcave].yIdx = id_map[primary_var[i]];

          //Now get gradient of i to get xIdx.
          int nnz;
          int & yIdx = simpleConcaves_[numberSimpleConcave].yIdx;
          int & xIdx = simpleConcaves_[numberSimpleConcave].xIdx;
          eval_grad_gi(n_var, NULL, false, i, nnz, jCol, NULL);
          if (nnz != 2) {//Error in ampl model
            std::cout<<"Incorrect suffixes description in ampl model. Constraint with id "
            <<id<<" is simple concave and should have only two nonzero elements"<<std::endl;
            exit(ERROR_IN_AMPL_SUFFIXES);
          }
          if (jCol[0] - 1== yIdx) {
            xIdx = jCol[1] - 1;
          }
          else {
            if (jCol[1] - 1!= yIdx) {//Error in ampl model
              std::cout<<"Incorrect suffixes description in ampl model. Constraint with id "
              <<id<<" : variable marked as y does not appear in the constraint."<<std::endl;
              exit(ERROR_IN_AMPL_SUFFIXES);
            }
            xIdx = jCol[0] - 1;
          }
          numberSimpleConcave++;
        }
      }
      delete [] jCol;
      numberSimpleConcave_ = numberSimpleConcave;
      numberNonConvex_ = numberSimpleConcave;
    }

  }


  void AmplTMINLP::read_onoff()
  {
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();
    DBG_ASSERT(asl);

    const AmplSuffixHandler * suffix_handler = GetRawPtr(suffix_handler_);
    const Index * onoff_c = suffix_handler->GetIntegerSuffixValues("onoff_c", AmplSuffixHandler::Constraint_Source);
    const Index * onoff_v = suffix_handler->GetIntegerSuffixValues("onoff_v", AmplSuffixHandler::Variable_Source);

    if(onoff_c == NULL && onoff_v == NULL){//No suffixes
      return;
    } 
    if(onoff_c == NULL || onoff_v == NULL){// If one in non-null both should be
        std::cerr<<"Incorrect suffixes description in ampl model.  One of per_v or per_c is declared but not the other."<<std::endl;
        exit(ERROR_IN_AMPL_SUFFIXES);
    } 

    c_extra_id_.clear();
    c_extra_id_.resize(n_con, -1);
    std::map<int, int> id_map;

      for (int i = 0 ; i < n_var ; i++) {
        if(onoff_v[i] > 0)
          id_map[onoff_v[i]] = i;
      }

      for (int i = 0 ; i < n_con ; i++) {
        if(onoff_c[i] > 0){
          std::map<int, int >::iterator k = id_map.find(onoff_c[i]);
          if(k != id_map.end()){
            c_extra_id_[i] = (*k).second;
          }
          else{
            std::cerr<<"Incorrect suffixes description in ampl model. onoff_c has value attributed to no variable "<<std::endl;
            exit(ERROR_IN_AMPL_SUFFIXES);
          }
        }
      }
   }

  bool AmplTMINLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
  {
    return ampl_tnlp_->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
  }

  bool AmplTMINLP::get_variables_types(Index n, VariableType* var_types)
  {
    // The variables are sorted by type in AMPL, so all we need to
    // know are the counts of each type.


    Index n_non_linear_b= 0;
    Index n_non_linear_bi= 0;
    Index n_non_linear_c= 0;
    Index n_non_linear_ci = 0;
    Index n_non_linear_o= 0;
    Index n_non_linear_oi = 0;
    Index n_binaries = 0;
    Index n_integers = 0;
    ampl_tnlp_->get_discrete_info(n_non_linear_b, n_non_linear_bi, n_non_linear_c,
        n_non_linear_ci, n_non_linear_o, n_non_linear_oi,
        n_binaries, n_integers);
    int totalNumberOfNonContinuous = 0;


    //begin with variables non-linear both in the objective function and in some constraints
    // The first ones are continuous
    Index start = 0;
    Index end = n_non_linear_b - n_non_linear_bi;
    for (int i=start; i<end; i++) {
      var_types[i] = CONTINUOUS;
    }

    // The second ones are integers
    start = end;
    end = start + n_non_linear_bi;
    for (int i=start; i < end; i++) {
      var_types[i] = INTEGER;
      totalNumberOfNonContinuous++;
    }
    //next variables non-linear in some constraints
    // The first ones are continuous
    start = end;
    end = std::max(start,end + n_non_linear_c - n_non_linear_ci - n_non_linear_b);
    for (int i=start; i<end; i++) {
      var_types[i] = CONTINUOUS;
    }

    // The second ones are integers
    start = end;
    end = start + n_non_linear_ci;
    for (int i=start; i < end; i++) {
      var_types[i] = INTEGER;
      totalNumberOfNonContinuous++;
    }

    //next variables non-linear in the objective function
    // The first ones are continuous
    start = end;
    end = std::max(start,end + n_non_linear_o - std::max(n_non_linear_b, n_non_linear_c) - n_non_linear_oi);
    for (int i=start; i<end; i++) {
      var_types[i] = CONTINUOUS;
    }

    // The second ones are integers
    start = end;
    end = start + n_non_linear_oi;
    for (int i=start; i < end; i++) {
      var_types[i] = INTEGER;
      totalNumberOfNonContinuous++;
    }

    //At last the linear variables
    // The first ones are continuous
    start = end;
    end = n - n_binaries - n_integers;
    for (int i=start; i<end; i++) {
      var_types[i] = CONTINUOUS;
    }

    // The second ones are binaries
    start = end;
    end = start + n_binaries;
    for (int i=start; i < end; i++) {
      var_types[i] = BINARY;
      totalNumberOfNonContinuous++;
    }

    // The third ones are integers
    start = end;
    end = start + n_integers;
    for (int i=start; i < end; i++) {
      var_types[i] = INTEGER;
      totalNumberOfNonContinuous++;
    }
    return true;
  }

  bool AmplTMINLP::get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types)
  {
    // The variables are sorted by type in AMPL, so all we need to
    // know are the counts of each type.


    Index n_non_linear_b= 0;
    Index n_non_linear_bi= 0;
    Index n_non_linear_c= 0;
    Index n_non_linear_ci = 0;
    Index n_non_linear_o= 0;
    Index n_non_linear_oi = 0;
    Index n_binaries = 0;
    Index n_integers = 0;
    ampl_tnlp_->get_discrete_info(n_non_linear_b, n_non_linear_bi, n_non_linear_c,
        n_non_linear_ci, n_non_linear_o, n_non_linear_oi,
        n_binaries, n_integers);

    //Compute the number of non linear variables:
    int n_non_linear = std::max(n_non_linear_c, n_non_linear_o);//n_non_linear_c + n_non_linear_o - n_non_linear_b;

    //printf("n_non_linear_c %i n_non_linear_o %i n_non_linear_b %i\n", n_non_linear_c, n_non_linear_o, n_non_linear_b);

    int start = 0;
    int end = n_non_linear;
    for (int i=start; i<end; i++) {
      var_types[i] = Ipopt::TNLP::NON_LINEAR;
    }

    //At last the linear variables
    // The first ones are continuous
    start = end;
    end = n;
    for (int i=start; i<end; i++) {
      var_types[i] = Ipopt::TNLP::LINEAR;
    }
    return true;
  }

  /** Returns the constraint linearity.
  * array should be alocated with length at least n.*/
  bool
  AmplTMINLP::get_constraints_linearity(Index m,
      Ipopt::TNLP::LinearityType* const_types)
  {
    return ampl_tnlp_->get_constraints_linearity(m, const_types);
  }

  bool AmplTMINLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
      Index m, Number* g_l, Number* g_u)
  {
    return ampl_tnlp_->get_bounds_info(n, x_l, x_u, m, g_l, g_u);
  }

  bool AmplTMINLP::get_starting_point(Index n, bool init_x, Number* x,
      bool init_z, Number* z_L, Number* z_U,
      Index m, bool init_lambda, Number* lambda)
  {
    return ampl_tnlp_->get_starting_point(n, init_x, x,
        init_z, z_L, z_U,
        m, init_lambda, lambda);
  }

  bool AmplTMINLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
  {
    return ampl_tnlp_->eval_f(n, x, new_x, obj_value);
  }

  bool AmplTMINLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
  {
    return ampl_tnlp_->eval_grad_f(n, x, new_x, grad_f);
  }

  bool AmplTMINLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
  {
    return ampl_tnlp_->eval_g(n, x, new_x, m, g);
  }

  bool AmplTMINLP::eval_jac_g(Index n, const Number* x, bool new_x,
      Index m, Index nele_jac, Index* iRow,
      Index *jCol, Number* values)
  {
    return ampl_tnlp_->eval_jac_g(n, x, new_x,
        m, nele_jac, iRow, jCol,
        values);
  }

  bool AmplTMINLP::eval_h(Index n, const Number* x, bool new_x,
      Number obj_factor, Index m, const Number* lambda,
      bool new_lambda, Index nele_hess, Index* iRow,
      Index* jCol, Number* values)
  {
    return ampl_tnlp_->eval_h(n, x, new_x, obj_factor,
        m, lambda, new_lambda, nele_hess, iRow,
        jCol, values);
  }

  bool AmplTMINLP::eval_gi(Index n, const Number* x, bool new_x,
      Index i, Number& gi)
  {
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();

    // ignore new_x for now
    xunknown();

    fint nerror = 0;
    gi = conival(i, const_cast<real*>(x), &nerror);
    if (nerror!=0) {
      return false;
    }
    else {
      return true;
    }
  }

  bool AmplTMINLP::eval_grad_gi(Index n, const Number* x, bool new_x,
      Index i, Index& nele_grad_gi, Index* jCol,
      Number* values)
  {
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();

    if (jCol) {
      // Only compute the number of nonzeros and the indices
      DBG_ASSERT(!values);
      nele_grad_gi = 0;
      for (cgrad* cg=Cgrad[i]; cg; cg = cg->next) {
        jCol[nele_grad_gi++] = cg->varno + 1;
      }
      return true;
    }
    DBG_ASSERT(values);

    // ignore new_x for now
    xunknown();

    asl->i.congrd_mode = 1;
    fint nerror = 0;
    congrd(i, const_cast<real*>(x), values, &nerror);
    if (nerror!=0) {
      return false;
    }
    else {
      return true;
    }
  }

  void AmplTMINLP::finalize_solution(TMINLP::SolverReturn status,
      Index n, const Number* x, Number obj_value)
  {
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();
    std::string message;
    std::string status_str;
    if (status == TMINLP::SUCCESS) {
      status_str = "\t\"Finished\"";
      message = "\n" + appName_ +": Optimal";
      solve_result_num = 3;
    }
    else if (status == TMINLP::INFEASIBLE) {
      status_str = "\t\"Finished\"";
      message = "\n" + appName_ + ": Infeasible problem";
      solve_result_num = 220;
    }
    else if (status == TMINLP::CONTINUOUS_UNBOUNDED) {
      status_str = "\t\"Finished\"";
      message = "\n" + appName_ +" Continuous relaxation is unbounded.";
      solve_result_num = 300;
    }
    else if (status == TMINLP::LIMIT_EXCEEDED) {
      status_str = "\t\"Not finished\"";
      message = "\n" + appName_ + ": Optimization interrupted on limit.";
      if(x)
        solve_result_num = 421; /* Limit reached or user interrupt with integer feasible solution.*/
      else
        solve_result_num = 410; /* Limit reached or user interrupt without solution.*/
    }
    else if (status == TMINLP::USER_INTERRUPT) {
      status_str = "\t\"Not finished\"";
      message = "\n" + appName_ + ": Optimization interrupted by user.";
      if(x)
        solve_result_num = 422; /* Limit reached or user interrupt with integer feasible solution.*/
      else
        solve_result_num = 411; /* Limit reached or user interrupt without solution.*/
    }
    else if (status == TMINLP::MINLP_ERROR) {
      status_str = "\t\"Aborted\"";
      message = "\n" + appName_ + ": Error encountered in optimization.";
      solve_result_num = 500;
    }
    if (writeAmplSolFile_) {
      write_solution(message, x);
      std::cout<<"\n "<<status_str<<std::endl;
    }
    else {
      std::cout<<status_str<<message<<std::endl;
      std::string fName = appName_ + ".sol";
      std::ofstream of(fName.c_str());
      for (int i = 0 ; i < n ; i++) {
        of<<i<<"\t"<<x[i]<<std::endl;
      }
      of<<"-1\n";
    }
  }

  void AmplTMINLP::write_solution(const std::string & message, const Number *x_sol)
  {
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();

    DBG_ASSERT(asl);
    //    DBG_ASSERT(x_sol);

    // We need to copy the message into a non-const char array to make
    // it work with the AMPL C function.
    char* cmessage = new char[message.length()+1];
    strcpy(cmessage, message.c_str());

    // In order to avoid casting into non-const, we copy the data here...
    double* x_sol_copy = NULL;
    if (x_sol) {
      x_sol_copy = new double[n_var];
      for (int i=0; i<n_var; i++) {
        x_sol_copy[i] = x_sol[i];
      }
    }
    write_sol(cmessage, x_sol_copy, NULL, NULL);

    delete [] x_sol_copy;
    delete [] cmessage;

  }


  /** This methods gives the linear part of the objective function */
  void AmplTMINLP::getLinearPartOfObjective(double * obj)
  {
    Index n, m, nnz_jac_g, nnz_h_lag;
    TNLP::IndexStyleEnum index_style = TNLP::FORTRAN_STYLE;
    get_nlp_info( n, m, nnz_jac_g, nnz_h_lag, index_style);
    eval_grad_f(n, NULL, 0,obj);
    Index n_non_linear_b= 0;
    Index n_non_linear_bi= 0;
    Index n_non_linear_c= 0;
    Index n_non_linear_ci = 0;
    Index n_non_linear_o= 0;
    Index n_non_linear_oi = 0;
    Index n_binaries = 0;
    Index n_integers = 0;
    ampl_tnlp_->get_discrete_info(n_non_linear_b, n_non_linear_bi, n_non_linear_c,
        n_non_linear_ci, n_non_linear_o, n_non_linear_oi,
        n_binaries, n_integers);

    // Now get the coefficients of variables wich are linear in the objective
    // The first ones are not
    Index start = 0;
    Index end = n_non_linear_b;
    for (int i=start; i<end; i++) {
      obj[i] = 0.;
    }

    //next variables should be linear in the objective just skip them
    // (from current end to (end + n_non_linear_c - n_non_linear_ci - n_non_linear_b;)


    // Those are nonlinear in the objective
    start = end + n_non_linear_c;
    end = start + n_non_linear_o;
    for (int i=start; i < end; i++) {
      obj[i]=0.;
    }
    //The rest is linear keep the values of the gradient
  }



  /** This method to returns the value of an alternative objective function for
  upper bounding (if one has been declared by using the prefix UBObj).*/
  bool
  AmplTMINLP::eval_upper_bound_f(Index n, const Number* x,
      Number& obj_value)
  {
    ASL_pfgh* asl = ampl_tnlp_->AmplSolverObject();
    //xknown(x);    // This tells ampl to use a new x
    fint nerror = -1;
    obj_value = objval(upperBoundingObj_, const_cast<double *>(x), &nerror);
    if (nerror > 0) {
      jnlst_->Printf(J_ERROR, J_MAIN,
          "Error in evaluating upper bounding objecting");
      throw -1;
    }
    return nerror;

  }
  /** Return the ampl solver object (ASL*) */
  const ASL_pfgh*
  AmplTMINLP::AmplSolverObject() const
  {
    return ampl_tnlp_->AmplSolverObject();
  }

} // namespace Bonmin

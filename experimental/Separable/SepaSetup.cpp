// (C) Copyright CNRS and others 2010
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, Université de la Méditérannée
// Hassan Hijazi, Orange Labs
//
// Date : 05/22/2010

#include "BonminConfig.h"
#include "OsiClpSolverInterface.hpp"

#include "SepaSetup.hpp"
#include "SepaTMINLP2OsiLP.hpp"
#include "SepaHeuristicInnerApproximation.hpp"
#include "BonOuterDescription.hpp"

namespace Sepa
{
  SepaSetup::SepaSetup(const CoinMessageHandler * handler):BonminSetup(handler)
  {
  }

  SepaSetup::SepaSetup(const SepaSetup &other):BonminSetup(other)
  {}

  SepaSetup::SepaSetup(const SepaSetup &other,
                           Bonmin::OsiTMINLPInterface &nlp):
      BonminSetup(other, nlp)
  {
  }

  SepaSetup::SepaSetup(const SepaSetup &other,
                           Bonmin::OsiTMINLPInterface &nlp,
                           const std::string &prefix):
    BonminSetup(other, nlp, prefix)
  {
   Bonmin::Algorithm algo = getAlgorithm();
    if (algo == Bonmin::B_OA)
      initializeSepa();
  }

  void SepaSetup::registerAllOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
     Bonmin::BonminSetup::registerAllOptions(roptions);

     Sepa::HeuristicInnerApproximation::registerOptions(roptions);

     roptions->SetRegisteringCategory("Initial Approximations descriptions", Bonmin::RegisteredOptions::UndocumentedCategory);
     roptions->AddStringOption2("initial_outer_description",
                                "Do we add all Outer Approximation constraints defining the initial Outer Approximation "
                                "description of the MINLP. See the number_approximations_initial_outer option for fixing the "
                                "number of approximation points",
                                "yes", "no","Do not generate the description", "yes","Generate the description",
                                "");

     roptions->AddUpperBoundedIntegerOption("number_approximations_initial_outer",
                                            "Number of Outer Approximation points needed for generating the initial Outer Approximation description, maximum value = 500, default value = 50",
                                            500, 50, "");
  }

  /** Register all the Bonmin options.*/
  void
  SepaSetup::registerOptions()
  {
    registerAllOptions(roptions_);
  }

  /** Initialize, read options and create appropriate bonmin setup using initialized tminlp.*/
  void
  SepaSetup::initialize(Ipopt::SmartPtr<Bonmin::TMINLP> tminlp, bool createContinuousSolver /*= false*/)
  {

    int do_outer;
    int n_approx;
    options()->GetEnumValue("initial_outer_description", do_outer, prefix_.c_str());
    options()->GetIntegerValue("number_approximations_initial_outer",
                               n_approx, prefix_.c_str());
    SepaTMINLP2OsiLP* linearizer = new SepaTMINLP2OsiLP;
    linearizer_ = linearizer;
    if(do_outer)
      linearizer->set_num_approx(n_approx);

    Bonmin::BonminSetup::initialize(tminlp, createContinuousSolver);

    if (getAlgorithm() == Bonmin::B_OA)
      initializeSepa();
  }

  /** Initialize, read options and create appropriate bonmin setup using initialized tminlp.*/
  void
  SepaSetup::initialize(const Bonmin::OsiTMINLPInterface &nlpSi, bool createContinuousSolver /*= false*/)
  {
    int do_outer;
    int n_approx;
    options()->GetEnumValue("initial_outer_description", do_outer, prefix_.c_str());
    options()->GetIntegerValue("number_approximations_initial_outer",
                               n_approx, prefix_.c_str());
    SepaTMINLP2OsiLP* linearizer = new SepaTMINLP2OsiLP;
    linearizer_ = linearizer;
    if(do_outer)
      linearizer->set_num_approx(n_approx);
    
    BonminSetup::initialize(nlpSi, createContinuousSolver);
    if (getAlgorithm() == Bonmin::B_OA)
      initializeSepa();
  }

  void SepaSetup::initializeSepa()
  {


    int doOuter;
    int nbAp = 10;
    options()->GetEnumValue("initial_outer_description", doOuter, prefix_.c_str());
    options()->GetIntegerValue("number_approximations_initial_outer",
                               nbAp, prefix_.c_str());

#ifdef USE_OLD_FUNC
    if(doOuter)
      addOuterDescription(*nonlinearSolver(), *continuousSolver(), nonlinearSolver()->getColSolution(), nbAp, false);
#endif
  
    int doInner;
    
    options()->GetEnumValue("heuristic_inner_approximation", doInner, prefix_.c_str());
    if(doInner){
      Sepa::HeuristicInnerApproximation * inner = new Sepa::HeuristicInnerApproximation(this);
      HeuristicMethod h;
      h.heuristic = inner;
      h.id = "InnerApproximation";
      heuristics_.push_back(h);
    }

  }

}/* end namespace Bonmin*/


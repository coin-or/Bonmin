// (C) Copyright CNRS and others 2010
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Université de la Méditérannée
// Hassan Hijazi, Orange Labs
//
// Date : 05/22/2010

#include "BonminConfig.h"
#include "OsiClpSolverInterface.hpp"

#include "SepaSetup.hpp"
#include "BonHeuristicInnerApproximation.hpp"
#include "BonOuterDescription.hpp"
namespace Bonmin
{
  SepaSetup::SepaSetup(const CoinMessageHandler * handler):BonminSetup(handler)
  {}

  SepaSetup::SepaSetup(const SepaSetup &other):BonminSetup(other)
  {}

  SepaSetup::SepaSetup(const SepaSetup &other,
                           OsiTMINLPInterface &nlp):
      BonminSetup(other, nlp)
  {
  }

  SepaSetup::SepaSetup(const SepaSetup &other,
                           OsiTMINLPInterface &nlp,
                           const std::string &prefix):
    BonminSetup(other, nlp, prefix)
  {
   Algorithm algo = getAlgorithm();
    if (algo == B_OA)
      initializeSepa();
  }

  void SepaSetup::registerAllOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
     BonminSetup::registerAllOptions(roptions);

     HeuristicInnerApproximation::registerOptions(roptions);

        roptions->SetRegisteringCategory("Initial Approximations descriptions", RegisteredOptions::UndocumentedCategory);
	roptions->AddStringOption2("initial_outer_description",
		"Do we add all Outer Approximation constraints defining the initial Outer Approximation description of the MINLP. See the number_approximations_initial_outer option for fixing the number of approximation points",
		"yes",
		"no","Do not generate the description",
		"yes","Generate the description",
		"");
	roptions->AddUpperBoundedIntegerOption("number_approximations_initial_outer",
		"Number of Outer Approximation points needed for generating the initial Outer Approximation description, maximum value = 500, default value = 50",
		500,
		50,
		"");
  }

  /** Register all the Bonmin options.*/
  void
  SepaSetup::registerOptions()
  {
    registerAllOptions(roptions_);
  }

  /** Initialize, read options and create appropriate bonmin setup using initialized tminlp.*/
  void
  SepaSetup::initialize(Ipopt::SmartPtr<TMINLP> tminlp, bool createContinuousSolver /*= false*/)
  {
    BonminSetup::initialize(tminlp, createContinuousSolver);
    if (getAlgorithm() == B_OA)
      initializeSepa();
  }

  /** Initialize, read options and create appropriate bonmin setup using initialized tminlp.*/
  void
  SepaSetup::initialize(const OsiTMINLPInterface &nlpSi, bool createContinuousSolver /*= false*/)
  {
    BonminSetup::initialize(nlpSi, createContinuousSolver);
    if (getAlgorithm() == B_OA)
      initializeSepa();
  }

  void SepaSetup::initializeSepa()
  {

    int doOuter;
    int nbAp = 10;
    options()->GetEnumValue("initial_outer_description", doOuter, prefix_.c_str());
    options()->GetIntegerValue("number_approximations_initial_outer",
       		nbAp, prefix_.c_str());
    if(doOuter)
      addOuterDescription(*nonlinearSolver(), *continuousSolver(), nonlinearSolver()->getColSolution(), nbAp, false);
    int doInner;
    
    options()->GetEnumValue("heuristic_inner_approximation", doInner, prefix_.c_str());
    if(doInner){
      HeuristicInnerApproximation * inner = new HeuristicInnerApproximation(this);
      HeuristicMethod h;
      h.heuristic = inner;
      h.id = "InnerApproximation";
      heuristics_.push_back(h);
    }

  }

}/* end namespace Bonmin*/


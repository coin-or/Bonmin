// (C) Copyright CNRS and others 2010
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, Université de la Méditérannée
//
// Date : 05/22/2010

#include <climits>

#include "BonminConfig.h"
#include "BonStrongBranchingSolver.hpp"

#include "RobotSetup.hpp"
#include "BonNWayObject.hpp"
#include "BonNWayChoose.hpp"



namespace Bonmin
{
  RobotSetup::RobotSetup(const CoinMessageHandler * handler):BonminSetup(handler)
  {}

  RobotSetup::RobotSetup(const RobotSetup &other):BonminSetup(other)
  {}

  RobotSetup::RobotSetup(const RobotSetup &other,
                           OsiTMINLPInterface &nlp):
      BonminSetup(other, nlp)
  {
  }

  RobotSetup::RobotSetup(const RobotSetup &other,
                           OsiTMINLPInterface &nlp,
                           const std::string &prefix):
    BonminSetup(other, nlp, prefix)
  {
    initializeRobot();
  }

  void RobotSetup::registerAllOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions)
  {
     BonminSetup::registerAllOptions(roptions);
     BonNWayChoose::registerOptions(roptions);


    roptions->AddLowerBoundedIntegerOption("branch_on_frac_only",
        "Starting at given depth branch on the subset of fractional variables (and set the last branch that one of them is 1)",
        0,INT_MAX,"");

    roptions->AddStringOption2("do_a_quick_one",
        "Do we try our luck?",
        "no",
        "no", "Don't (of course).",
        "yes", "Be crazy",
        "");

  }

  /** Register all the Bonmin options.*/
  void
  RobotSetup::registerOptions()
  {
    registerAllOptions(roptions_);
  }

  /** Initialize, read options and create appropriate bonmin setup using initialized tminlp.*/
  void
  RobotSetup::initialize(Ipopt::SmartPtr<TMINLP> tminlp, bool createContinuousSolver /*= false*/)
  {
    BonminSetup::initialize(tminlp,createContinuousSolver);
    initializeRobot();
  }

  /** Initialize, read options and create appropriate bonmin setup using initialized tminlp.*/
  void
  RobotSetup::initialize(const OsiTMINLPInterface &nlpSi, bool createContinuousSolver /*= false*/)
  {
    BonminSetup::initialize(nlpSi,createContinuousSolver);
    initializeRobot();
  }

  void RobotSetup::initializeRobot()
  {
    assert(continuousSolver_ == nonlinearSolver_);

    delete branchingMethod_;

    continuousSolver_->deleteObjects();
    continuousSolver_->findIntegersAndSOS(false);
    setPriorities();
    addNWays();
    Ipopt::SmartPtr<StrongBranchingSolver> strong_solver = NULL;
    nonlinearSolver_->SetStrongBrachingSolver(strong_solver);
    BonNWayChoose* chooseVariable = new BonNWayChoose(*this,
                                                      nonlinearSolver_);
    branchingMethod_ = chooseVariable;
    branchingMethod_->setNumberStrong(intParam_[NumberStrong]);

  }

  void
  RobotSetup::addNWays()
  {

    int do_quick;
    options()->GetEnumValue("do_a_quick_one", do_quick, prefix());
    int depth_frac;
    options()->GetIntegerValue("branch_on_frac_only", depth_frac, prefix());

    // pass user set Sos constraints (code inspired from CoinSolve.cpp)
    const TMINLP::SosInfo * sos = nonlinearSolver()->model()->sosConstraints();
    if (!getIntParameter(BabSetupBase::DisableSos) && sos && sos->num > 0) //we have some sos constraints
    {
      const int & numSos = sos->num;
      OsiObject ** objects = new OsiObject*[numSos];
      const int * starts = sos->starts;
      const int * indices = sos->indices;
      //const char * types = sos->types;
      const double * weights = sos->weights;
      bool hasPriorities = false;
      const int * varPriorities = nonlinearSolver()->getPriorities();
      int numberObjects =  nonlinearSolver()->numberObjects();
      if (varPriorities)
      {
        for (int i = 0 ; i < numberObjects ; i++) {
          if (varPriorities[i]) {
            hasPriorities = true;
            break;
          }
        }
      }
      const int * sosPriorities = sos->priorities;
      if (sosPriorities)
      {
        for (int i = 0 ; i < numSos ; i++) {
          if (sosPriorities[i]) {
            hasPriorities = true;
            break;
          }
        }
      }

      std::vector<std::list<int> > groups(numSos + 1);

      for (int i = 0 ; i < numSos ; i++)
      {
        int start = starts[i];
        int length = starts[i + 1] - start;
          for(int j = 0 ; j < length ; j++){
              groups[(size_t) weights[j]].push_back(indices[start+j]);
          }
      }

      for (int i = 0 ; i < numSos ; i++)
      {
        int start = starts[i];
        int length = starts[i + 1] - start;
          BonNWayObject * nway = new BonNWayObject(length, &indices[start],i);
          nway->setPriority(1);
          for(int j = 0 ; j < length ; j++){//Setup consequences
             n_way_consequences cons;
             std::vector<int>& ids = cons.indices;
             int idx = (int) weights[j];
             const std::list<int> &to_add = groups[idx];
             for(std::list<int>::const_iterator k = to_add.begin() ; 
                 k != to_add.end() ; k++){
               if(*k != indices[start+j]) ids.push_back(*k);
             }
           nway->setConsequence(j, cons);
          }
          objects[i] = nway;

        if(do_quick)
          nway->make_quick();
        nway->set_only_frac_branches(depth_frac);
        if (hasPriorities && sosPriorities && sosPriorities[i]) {
          objects[i]->setPriority(sosPriorities[i]);
        }
      }
      nonlinearSolver()->addObjects(numSos, objects);
      for (int i = 0 ; i < numSos ; i++)
        delete objects[i];
      delete [] objects;
    }
  }
}/* end namespace Bonmin*/


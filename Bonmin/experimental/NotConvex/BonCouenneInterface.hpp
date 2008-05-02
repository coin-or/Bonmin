// (C) Copyright International Business Machines Corporation (IBM) 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pietro Belotti, Carnegie Mellon University
// Pierre Bonami, International Business Machines Corporation
//
// Date : 12/19/2006

#ifndef CouenneInterface_H
#define CouenneInterface_H
#include "BonAmplInterface.hpp"
#include "CouenneCutGenerator.hpp"


struct ASL;

struct ASL *readASLfg (char **);


namespace Bonmin {
class CouenneInterface : public AmplInterface
{
  public:
  /** Default constructor. */
  CouenneInterface();

  /** Copy constructor. */
  CouenneInterface(const CouenneInterface &other);

  /** virutal copy constructor. */
  virtual CouenneInterface * clone(bool CopyData);

  /** Destructor. */
  virtual ~CouenneInterface();

    
  /** read ampl file using inputed options, journalist....*/
  virtual void readAmplNlFile(char **& argv, Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions,
                              Ipopt::SmartPtr<Ipopt::OptionsList> options,
                              Ipopt::SmartPtr<Ipopt::Journalist> journalist);
  
  /** \name Overloaded methods to build outer approximations */
  //@{
  /** \brief Extract a linear relaxation of the MINLP.
   * Solve the continuous relaxation and takes first-order 
   * outer-approximation constraints at the optimum.
   * The put everything in an OsiSolverInterface.
   */
  virtual void extractLinearRelaxation
  (OsiSolverInterface &si,  CouenneCutGenerator & couenneCg, bool getObj = 1, bool solveNlp = 1);

  
  /** To set some application specific defaults. */
  virtual void setAppDefaultOptions(Ipopt::SmartPtr<Ipopt::OptionsList> Options);

  /// return value of have_nlp_solution_
  bool haveNlpSolution ()
  {return have_nlp_solution_;}

protected:

  /// true if we got an integer feasible solution from initial solve 
  bool have_nlp_solution_;
};

} /** end Bonmin namespace. */
#endif

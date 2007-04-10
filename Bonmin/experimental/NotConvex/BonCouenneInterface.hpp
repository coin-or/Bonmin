// (C) Copyright International Business Machines Corporation (IBM) 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pietro Belloti, Carnegie Mellon University
// Pierre Bonami, International Business Machines Corporation
//
// Date : 12/19/2006

#ifndef CouenneInterface_H
#define CouenneInterface_H
#include "BonAmplInterface.hpp"
#include "CouenneCutGenerator.h"


struct ASL;

struct ASL *readASLfg (char **);


namespace Bonmin {
class CouenneInterface : public AmplInterface
{
  public:
  /** Default constructor. */
  CouenneInterface();

  /** Constructor with inputed ampl command line.*/
  CouenneInterface(char **& amplArgs, SmartPtr<TNLPSolver> app);

  /** Copy constructor. */
  CouenneInterface(const CouenneInterface &other);

  /** virutal copy constructor. */
  virtual CouenneInterface * clone(bool CopyData);

  /** Destructor. */
  virtual ~CouenneInterface();

  /** \name Overloaded methods to build outer approximations */
  //@{
  /** \brief Extract a linear relaxation of the MINLP.
   * Solve the continuous relaxation and takes first-order 
   * outer-approximation constraints at the optimum.
   * The put everything in an OsiSolverInterface.
   */
  virtual void extractLinearRelaxation(OsiSolverInterface &si, bool getObj = 1, bool solveNlp = 1);

  /** Get the outer approximation constraints at the currently stored optimal point.
   (Only get outer-approximations of nonlinear constraints of the problem.)*/
  //  virtual void getOuterApproximation(OsiCuts &cs, bool getObj, const double * x2, bool global);

  /** Get the outer approximation constraints at provided point.
      If x2 is different from NULL only add cuts violated by x2.
   (Only get outer-approximations of nonlinear constraints of the problem.)*/
  //  virtual void getOuterApproximation(OsiCuts &cs, const double * x, bool getObj, const double * x2, bool global);

  // This is never called, eliminated
  //  void updateCouenneAuxiliryVar(double *x, double * lb, double * ub);

  const CouenneProblem * couenneProb() const{
    return couenneCg_->Problem();
  }

  CouenneProblem * couenneProb() {
    return couenneCg_->Problem();
  }
  
  /** return ASL interface */
  ASL *getASL() {return aslfg_;}

  CouenneCutGenerator * couenneCg()  {
    return couenneCg_;}

  const CouenneCutGenerator * couenneCg()  const{
    return couenneCg_;}

  /** To set some application specific defaults. */
  virtual void setAppDefaultOptions(Ipopt::SmartPtr<Ipopt::OptionsList> Options);

protected:
 
  /** The cut generator from couenne. */
  CouenneCutGenerator *couenneCg_;

  /** The simpler ASL interface (no partially separable information,
      no groups, no hassle) */
  ASL *aslfg_;
};

} /** end Bonmin namespace. */
#endif

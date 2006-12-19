// (C) Copyright International Business Machines Corporation (IBM and Carnegie Mellon University 2006 
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 02/19/2006
#ifndef CouenneInterface_H
#define CouenneInterface_H
#include "BonAmplInterface.hpp"
#include "CouenneCutGenerator.h"


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
  virtual OsiSolverInterface * clone(bool CopyData);

  /** Destructor. */
  virtual ~CouenneInterface();

  /** \name Overloaded methods to build outer approximations */
  //@{
  /** \brief Extract a linear relaxation of the MINLP.
   * Solve the continuous relaxation and takes first-order outer-approximation constraints at the optimum.
   * The put everything in an OsiSolverInterface.
   */
  virtual void extractLinearRelaxation(OsiSolverInterface &si, bool getObj = 1);

  /** Get the outer approximation constraints at the currently stored optimal point.
   (Only get outer-approximations of nonlinear constraints of the problem.)*/
  virtual void getOuterApproximation(OsiCuts &cs, bool getObj = 1);

  protected:
  /** The cut generator from couenne. */
  CouenneCutGenerator *couenneCg_;

};


} /** end Bonmin namespace. */
#endif


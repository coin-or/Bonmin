 // (C) Copyright International Business Machines (IBM) 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, International Business Machines
//
// Date :  12/20/2006

#ifndef BonCouenneConvexCuts_HPP
#define BonCouenneConvexCuts_HPP

#include "BonOaDecBase.hpp"
#include "CglCutGenerator.hpp"
namespace Bonmin {
  class CouenneConvCuts: public OaDecompositionBase {
    public:
    CouenneConvCuts(OsiTMINLPInterface *nlp = NULL,
            int numRounds = 1
	    ):
      OaDecompositionBase(nlp,NULL, NULL,0,0,0),
      objValue_(-DBL_MAX),
      numRounds_(numRounds){
    }
    
    /// Copy constructor
    CouenneConvCuts(const CouenneConvCuts & copy):
      OaDecompositionBase(copy),
      objValue_(copy.objValue_),
      numRounds_(copy.numRounds_)
      {
    }
    
    /// clone
    CglCutGenerator * clone() const{
      return new CouenneConvCuts(*this);
    }
 
    /// Destructor
    virtual ~CouenneConvCuts(){
    }
    /** Standard cut generation methods. */
    virtual void generateCuts(const OsiSolverInterface &si,  OsiCuts & cs,
                              const CglTreeInfo info = CglTreeInfo()) const;
    double doCouenneConvRounds(OsiSolverInterface &si,
                     bool leaveSiUnchanged);

   void setNumRounds(int value){
   numRounds_ = value;}
   protected: 
    /// virtual method which performs the OA algorithm by modifying lp and nlp.
    virtual double performOa(OsiCuts & cs, solverManip &nlpManip, solverManip &lpManip, 
                           SubMipSolver *& subMip, OsiBabSolver * babInfo, double &cutoff) const{
    throw -1;
    }
    /// virutal method to decide if local search is performed
    virtual bool doLocalSearch() const{
      return 0;}
    private:
      /** Record obj value at final point of CouenneConv. */
      mutable double objValue_;
     /** number of iterations of generation. */
     int numRounds_;
  };
} /* end namespace Bonmin.*/
#endif

 // (C) Copyright International Business Machines (IBM) 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, International Business Machines
//
// Date :  12/20/2006

#ifndef BonECPCuts_HPP
#define BonECPCuts_HPP

#include "BonOaDecHelper.hpp"
#include "CglCutGenerator.hpp"
namespace Bonmin {
  class EcpCuts: public OaDecompositionHelper, public CglCutGenerator {
    public:
    EcpCuts(OsiTMINLPInterface *nlp = NULL
	    ):
      OaDecompositionHelper(nlp,NULL, NULL,0,0,0){
    }
    
    /// Copy constructor
    EcpCuts(const EcpCuts & copy):
      OaDecompositionHelper(copy),
      CglCutGenerator(copy){
    }
    
    /// clone
    CglCutGenerator * clone() const{
      return new EcpCuts(*this);
    }
 
    /// Destructor
    virtual ~EcpCuts(){
    }
    
    virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
			       const CglTreeInfo info = CglTreeInfo()) const;
    
  };
} /* end namespace Bonmin.*/
#endif

// (C) Copyright CNRS 2009
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, LIF, CNRS 
//
// Date : 13/02/2009

#ifndef BonCutStorage_H
#define BonCutStorage_H

#include "BonminConfig.h"
#include "CglCutGenerator.hpp"

namespace Bonmin {
  /** Class used to store cuts found in a subalgorithm which may be added in next cut
      generation round. */
  class CutStorage : public CglCutGenerator {
    public:
      /** Standard cut generation method. */
      void generateCuts(const OsiSolverInterface & si, OsiCuts &csi,
        const CglTreeInfo info = CglTreeInfo()) const;

      /** Add some cuts to the collection.*/
      void addSomeCuts(const OsiCut ** cuts);

    private:
    /** Stores the current cut*/
    mutable OsiCuts cs_; 
  } ;

} /** Ends Namespace Bonmin.*/

#endif


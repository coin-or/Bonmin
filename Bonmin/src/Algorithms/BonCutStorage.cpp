// (C) Copyright CNRS 2009
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, LIF, CNRS 
//
// Date : 13/02/2009

#include "BonCutStorage.hpp"

namespace Bonmin {

void
CutStorage::generateCuts(const OsiSolverInterface & si, OsiCuts &csi,
                         const CglTreeInfo info) const{
   int numRowCuts = cs_.sizeRowCuts();
   for(int i = 0 ; i != numRowCuts ; i++){
     OsiRowCut * cut = cs_.rowCutPtr(i);
     csi.insert(cut);
   }
   cs_.dumpCuts();
}

} /** Ends Namespace Bonmin.*/

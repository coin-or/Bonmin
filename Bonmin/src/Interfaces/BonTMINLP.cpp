// (C) Copyright International Business Machines (IBM) 2005, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, IBM
//
// Date : 26/09/2006

#include "BonTMINLP.hpp"
#include "IpBlas.hpp"

namespace Bonmin{

/** default constructor for Sos constraints */
TMINLP::SosInfo::SosInfo():
        num(0), 
        types(NULL), 
        priorities(NULL), 
        numNz(0), 
        starts(NULL),
        indices(NULL), 
        weights(NULL)
{}

/** Copy constructor.*/
TMINLP::SosInfo::SosInfo(const SosInfo & source):
        num(source.num), 
        types(NULL), 
        priorities(NULL), 
        numNz(source.numNz), 
        starts(NULL),
        indices(NULL),
        weights(NULL)
{

  if(num > 0) {
    assert(source.types!=NULL);
    assert(source.priorities!=NULL);
    assert(source.starts!=NULL);
    assert(source.indices!=NULL);
    assert(source.weights!=NULL);
    types = new char[num];
    priorities = new int[num];
    starts = new int[num + 1];
    indices = new int[numNz];
    weights = new double[numNz];
    for(int i = 0 ; i < num ; i++) {
      source.types[i] = types[i];
      source.priorities[i] = priorities[i];
      source.starts[i] = starts[i];
    }
    for(int i = 0 ; i < numNz ; i++) {
      source.indices[i] = indices[i];
      source.weights[i] = weights[i];
    }
  }
  else {
    assert(source.types==NULL);
    assert(source.priorities==NULL);
    assert(source.starts==NULL);
    assert(source.indices==NULL);
    assert(source.weights==NULL);
  }

}


/** Reset information */
void 
TMINLP::SosInfo::gutsOfDestructor()
{
  num = 0;
  numNz = 0;
  if(types) delete [] types;
  types = NULL;
  if(starts) delete [] starts;
  starts = NULL;
  if(indices) delete [] indices;
  indices = NULL;
  if(priorities) delete [] priorities;
  priorities = NULL;
  if(weights) delete [] weights;
  weights = NULL;
}


void TMINLP::PerturbInfo::SetPerturbationArray(Index numvars, const double* perturb_radius) {
  delete [] perturb_radius_;
  if (perturb_radius) {
    perturb_radius_ = new double[numvars];
    for(int i=0; i<numvars; i++) {
      perturb_radius_[i] = perturb_radius[i];
    }
  }
}

TMINLP::TMINLP()
{}

TMINLP::TMINLP(const TMINLP & source)
{
}
     
/** Default destructor */
TMINLP::~TMINLP()
{
}


}

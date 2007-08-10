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

TMINLP::TMINLP():
    cutsjCol_(NULL),
    cutsiRow_(NULL),
    cutsElems_(NULL),
    cutsLower_(NULL),
    cutsUpper_(NULL),
    nLinearCuts_(0),
    linearCutsNnz_(0),
    linearCutsCapacity_(0),
    linearCutsNnzCapacity_(0) 
    {}


/** Default destructor */
virtual 
TMINLP::~TMINLP()
{
   delete [] cutsElems_;
   delete [] cutsiRow_;
   delete [] cutsjCol_;
   delete [] cutsLower_;
   delete [] cutsUpper_;
}


void
TMINLP::addCuts(int numberCuts, const OsiRowCut ** cuts){

 int n,m,nnz_lag,nnz_hess;
 Ipopt::TNLP::IndexStyleEnum fort;
 get_nlp_info(n,m,nnz_lag,nnz_hess,fort);

 //count the number of non-zeroes
 //Cuts have to be added by rows for removeCuts to work
 int nnz = linearCutsNnz_;
 for(int i = 0 ; i < numberCuts ; i++)
 {
   nnz += cuts[i]->row().getNumElements();
 }
 int resizeNnz = nnz > linearCutsNnzCapacity_;
 int resizeCuts = nLinearCuts_ + numberCuts > linearCutsCapacity_;
 resizeLinearCuts(max((1+ resizeCuts) * linearCutsCapacity_, resizeCuts * numberCuts + linearCutsCapacity_), 
                  max((1 + resizeNnz) * linearCutsNnzCapacity_, resizeNnz * nnz + linearCutsNnzCapacity_));

 /* reinit nnz */
 nnz = linearCutsNnz_; 

 for(int i = 0 ; i < numberCuts ; i++)
 {
   const int * ind = cuts[i]->row().getIndices();
   const double * values = cuts[i]->row().getElements();
   int size = cuts[i]->row().getNumElements();
   int iplusn = i + nLinearCuts_;
   for(int j = 0; j < size ; j++)
   {
    DBG_ASSERT(ind[j] < n);
    cutsjCol_[nnz] = ind[j];
    cutsiRow_[nnz] = iplusn;
    cutsElems_[nnz++] = values[j];
   }
   cutsLower_[iplusn] = cuts[i]->lb();
   cutsUpper_[iplusn] = cuts[i]->ub();
  }
 nLinearCuts_+=numberCuts;
 linearCutsNnz_=nnz;
}

void
TMINLP::resizeLinearCuts(int newNumberCuts, int newNnz)
{
  if(newNumberCuts > linearCutsCapacity_)
  {
     double * newLower = new double[newNumberCuts];
     double * newUpper = new double[newNumberCuts];
     if(linearCutsCapacity_)
     {
       IpBlasDcopy(nLinearCuts_, cutsLower_, 1, newLower, 1);
       IpBlasDcopy(nLinearCuts_, cutsUpper_, 1, newUpper, 1);
       delete [] cutsLower_;
       delete [] cutsUpper_;
     }
     cutsLower_ = newLower;
     cutsUpper_ = newUpper;
     linearCutsCapacity_ = newNumberCuts;
  }
  if(newNnz > linearCutsNnzCapacity_)
  {
    double * newElems = new double [newNnz];
    int * newiRow = new int [newNnz];
    int * newjCol = new int [newNnz];
    if(linearCutsNnzCapacity_)
    {
      IpBlasDcopy(linearCutsNnz_, cutsElems_, 1, newElems, 1);
      for(int i = 0 ; i < linearCutsNnz_ ; i++) 
      {
        newiRow[i] = cutsiRow_[i];
        newjCol[i] = cutsjCol_[i];
      }
      delete [] cutsElems_;
      delete [] cutsiRow_;
      delete [] cutsjCol_;
    }
    cutsElems_ = newElems;
    cutsiRow_ = newiRow;
    cutsjCol_ = newjCol;
    linearCutsNnzCapacity_ = newNnz;
  }
}
void
TMINLP::removeCuts(int number, const int * toRemove){
   if(number==0) return;
    int n,m,nnz_lag,nnz_hess;
   Ipopt::TNLP::IndexStyleEnum fort;
   get_nlp_info(n,m,nnz_lag,nnz_hess,fort);

   int * sorted = new int[number];
   for(int i = 0 ; i < number ; i++) sorted[i] = toRemove[i] - m;
   std::sort(sorted, sorted + number);  
   int iNew = 0;
   int k = 0;
   for(int i = 0 ; i < linearCutsNnz_ ; i++)
   {
     if(sorted[k] < cutsiRow_[i]) k++;
     if(sorted[k] < cutsiRow_[i]) throw -1;
     if(sorted[k] == cutsiRow_[i]) continue;
     cutsiRow_[iNew] = cutsiRow_[i] - k;
     cutsjCol_[iNew] = cutsjCol_[i];
     cutsElems_[iNew++] = cutsElems_[i];
   }
 k=0;
 for(int i = sorted[k] ; i < nLinearCuts_ ; i++){
    if(sorted[k] < i) k++;
    if(sorted[k] < i) throw -1;
    if(sorted[k] == i) continue;
    cutsLower_[i - k] = cutsLower_[i];
    cutsUpper_[i - k] = cutsUpper_[i];
 }
 linearCutsNnz_ = iNew;
 nLinearCuts_ -= number;
 delete [] sorted;
}
void
TMINLP::removeLastCuts(int number){
 number = nLinearCuts_ - number;
 int iNew = 0;
   for(int i = 0 ; i < linearCutsNnz_ ; i++)
   {
     if(cutsiRow_[i] >= number) {
     continue;}
     cutsiRow_[iNew] = cutsiRow_[i];
     cutsjCol_[iNew] = cutsjCol_[i];
     cutsElems_[iNew++] = cutsElems_[i];
   }
 linearCutsNnz_ = iNew;
 nLinearCuts_ = number;
 std::cout<<"Number of cuts remaining "<<nLinearCuts_<<std::endl;
 }
}

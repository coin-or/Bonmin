// (C) Copyright International Business Machines Corporation, 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 11/21/2006




#ifndef BonFilterWarmStart_H
#define BonFilterWarmStart_H

#include "CoinWarmStartBasis.hpp"
#include "BonFilterSolver.hpp" /* for types */

#include <vector>
namespace Bonmin{

  /** Warm start for filter interface.
  Warm start for filter constists of a (possibly huge) array of integers.
  \bug Inheritance from CoinWarmStartBasis is only for compatibility with Cbc
  */
  class FilterWarmStart : public CoinWarmStartBasis
  {
    typedef FilterSolver::fint fint;
    typedef FilterSolver::real real;
    
  public:
    /** Default values for istat */
    static fint def_istat[14];
    /** Constructor */
    FilterWarmStart(const fint size = 0, const fint * warmArray = NULL, const fint istat[14] = def_istat):
      CoinWarmStartBasis(),
      size_(size),
      warmArray_(NULL){
      if (size_ > 0){
	warmArray_ = new fint[size];
	if(warmArray != NULL){
	  CoinCopyN(warmArray, size, warmArray_);
	}
      }
      else if(warmArray != NULL) {
	throw CoinError("Array passed but size is 0","FilterWarmStart(const fint, const fint *)","FilterWarmStart");
	}
      for(int i = 0 ; i < 14 ; i ++)
	istat_[i] = istat[i];
    }

    /** Copy constructor */
    FilterWarmStart(const FilterWarmStart & other):
      CoinWarmStartBasis(other),
      size_(other.size_),
      warmArray_(NULL){
      if (size_ > 0){
	warmArray_ = new fint[size_];
	if(other.warmArray_ != NULL){
	  CoinCopyN(other.warmArray_, size_, warmArray_);
	}
      }
      else if(other.warmArray_ != NULL) {
	throw CoinError("Array passed but size is 0","FilterWarmStart(const fint, const fint *)","FilterWarmStart");
	}      
    }

    /** virtual copy */
    virtual CoinWarmStart * clone() const
    { return new FilterWarmStart(*this);}

    /**Set size of the array. */
    void setInfo(const fint size = 0, const fint * warmArray = NULL, const fint istat[14] = def_istat){
      if(size != size_){
	size_ = size;
	if(warmArray_) delete [] warmArray_;
	warmArray_ = NULL;
	if(size > 0)
	  warmArray_ = new fint[size];
      }
      else if(size > 0){
	assert(warmArray_);
      }
      if(size <= 0 && warmArray)
	throw CoinError("Array passed but size is 0","setInfo(const fint, const fint *)","FilterWarmStart");	  
      CoinCopyN(warmArray, size_, warmArray_);
      
      for(int i = 0 ; i < 14 ; i ++)
	istat_[i] = istat[i];
    }
    
    /** Destructor. */
    virtual ~FilterWarmStart(){
      delete [] warmArray_;
    }

    /** Generate differences.*/
    virtual CoinWarmStartDiff * generateDiff(const CoinWarmStart * const other) const;

    /** Apply differences. */
    virtual void applyDiff(const CoinWarmStartDiff * const cswDiff);


    /** Access to array */
    const fint * array() const{
      return warmArray_;
    }

    /** Access to size. */
    fint size() const {
      return size_;}

    const fint * istat()const {
      return istat_;}
  private:
    /** Size of the warm start object. */
    fint size_;

    /** Warm start information */
    fint * warmArray_;

    /** istat */
    fint istat_[14];
    
  };


  class FilterWarmStartDiff : public CoinWarmStartBasisDiff
  {
    typedef FilterSolver::fint fint;
    typedef FilterSolver::real real;

  public:
    FilterWarmStartDiff(fint capacity):
      CoinWarmStartBasisDiff(),
      differences(0){
      differences.reserve(capacity);
    }

    virtual CoinWarmStartDiff * clone() const{
      int size = differences.size();
      FilterWarmStartDiff * return_value = new FilterWarmStartDiff(size);
      return_value->differences = differences;
      for(int i = 0 ; i < 14 ; i++)
	{
	  return_value->istat_[i] = istat_[i];
	}
      return return_value;
    }
    friend class FilterWarmStart;
  private:
    /** One difference is two integers (indice and difference). */
    typedef std::pair<fint, fint> OneDiff;
    /** Vector of all the differences.*/
    std::vector<OneDiff> differences;

    /** istat */
    fint istat_[14];

  };

} /* end namespace Bonmin */
#endif


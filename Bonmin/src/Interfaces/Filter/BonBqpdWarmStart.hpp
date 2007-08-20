// (C) Copyright International Business Machines Corporation 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Andreas Waechter, International Business Machines Corporation
//                    based on BonFilterWarmStart.hpp
//
// Date : 2007-08-03


#ifndef BonBqpdWarmStart_H
#define BonBqpdWarmStart_H

#include "CoinWarmStartBasis.hpp"
#include "BonFilterSolver.hpp" /* for types */

#include <vector>

namespace Bonmin{

  /** Warm start for filter interface.
  Warm start for filter constists of a (possibly huge) array of integers.
  \bug Inheritance from CoinWarmStartBasis is only for compatibility with Cbc
  */
  class BqpdWarmStart : public CoinWarmStartBasis
  {
    typedef FilterSolver::fint fint;
    typedef FilterSolver::real real;
    
  public:
    /** Default values for istat */
    static fint def_istat[14];
    /** Constructor */
    BqpdWarmStart(const fint xSize = 0,
		    const real* xArray = NULL,
		    const fint lamSize = 0,
		    const real* lamArray = NULL,
		    const fint lwsSize = 0,
		    const fint *lwsArray = NULL,
		    const fint istat[14] = def_istat);

    /** Copy constructor */
    BqpdWarmStart(const BqpdWarmStart & other);

    /** virtual copy */
    virtual CoinWarmStart * clone() const
    { return new BqpdWarmStart(*this);}

#ifdef AWDoesntseemnecessary
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
	throw CoinError("Array passed but size is 0","setInfo(const fint, const fint *)","BqpdWarmStart");	  
      CoinCopyN(warmArray, size_, warmArray_);
      
      for(int i = 0 ; i < 14 ; i ++)
	istat_[i] = istat[i];
    }
#endif

    /** Destructor. */
    virtual ~BqpdWarmStart();

    /** Generate differences.*/
    virtual CoinWarmStartDiff* generateDiff(const CoinWarmStart * const other) const;

    /** Apply differences. */
    virtual void applyDiff(const CoinWarmStartDiff * const cswDiff);

    /** Access to x Array */
    const real *xArray() const {
      if (tempxArray_) {
	return tempxArray_;
      }
      else {
	return xArray_;
      }
    }

    /** Array to x size */
    fint xSize() const {
      return xSize_;}

    /** Access to lam Array */
    const real *lamArray() const {
      if (templamArray_)
	return templamArray_;
      else
	return lamArray_;
    }

    /** Array to lam size */
    fint lamSize() const {
      return lamSize_;}

    /** Access to lws array */
    const fint *lwsArray() const{
      return lwsArray_;
    }

    /** Access to lws size. */
    fint lwsSize() const {
      return lwsSize_;}

    const fint* istat()const {
      return istat_;}

    void flushPoint();
  private:
    /** Size of real x array store. */
    fint xSize_;

    /** Real x array to store */
    real* xArray_;

    /** Size of real lam array store. */
    fint lamSize_;

    /** Real lam array to store */
    real* lamArray_;

    /** Size of fint lws array store. */
    fint lwsSize_;

    /** fint lws array to store */
    fint* lwsArray_;

    /** Filter's istat (AW: I think we only need first entry) */
    fint istat_[14];
  };

} /* end namespace Bonmin */
#endif


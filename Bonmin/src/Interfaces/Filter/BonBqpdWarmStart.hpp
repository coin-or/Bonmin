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
#include "CoinWarmStartPrimalDual.hpp"
#include "BonFilterSolver.hpp" /* for types */

#include <vector>

namespace Bonmin{

  /** Warm start for filter interface.
   * Warm start for filter constists of a (possibly huge) array of integers.
   * This class inherits from CoinWarmStartPrimalDual, because that's what
   * this warmstart really is. <br>
   * For practical reason (integration in Cbc) this class also inherits from
   * CoinWarmStartBasis. <br>
   */
  class BqpdWarmStart :
    public virtual CoinWarmStartPrimalDual, public virtual CoinWarmStartBasis,
    public Coin::ReferencedObject
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

    /** constructor from a CoinWarmStartPrimalDual */
    BqpdWarmStart(const CoinWarmStartPrimalDual& pdws);

    /** virtual copy */
    virtual CoinWarmStart * clone() const
    { return new BqpdWarmStart(*this);}

    /** Destructor. */
    virtual ~BqpdWarmStart();

    /** Generate differences.*/
    virtual CoinWarmStartDiff* generateDiff(const CoinWarmStart * const other) const;

    /** Apply differences. */
    virtual void applyDiff(const CoinWarmStartDiff * const cswDiff);

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

    ///Is this an empty warm start?
    bool empty() const
    {
      return empty_;
    }
  private:
    /** Size of fint lws array store. */
    fint lwsSize_;

    /** fint lws array to store */
    fint* lwsArray_;

    /** Filter's istat (AW: I think we only need first entry) */
    fint istat_[14];
    ///Say if warm start is empty
    bool empty_;
  };

} /* end namespace Bonmin */
#endif


// (C) Copyright International Business Machines Corporation, Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Carnegie Mellon University,
// Andreas Waechter, International Business Machines Corporation
//
// Date : 02/15/2006


#ifndef IpoptWarmStart_HPP
#define IpoptWarmStart_HPP
#include "CoinWarmStartBasis.hpp"
#include "CoinPackedVector.hpp"
#include "BonIpoptInteriorWarmStarter.hpp"


namespace Bonmin {
class TMINLP2TNLP;

/** \brief Class for storing warm start informations for Ipopt.<br>
 * For practical reason (integration in Cbc) this class inherits from
 * CoinWarmStartBasis. <br>
 * This class stores a starting point (primal and dual values) for Ipopt.
 */
class IpoptWarmStart : public CoinWarmStartBasis
{
public:

  /// Default constructor
  IpoptWarmStart(bool empty = 1, int numvars = 0, int numcont = 0);
  /// Usefull constructor, stores the current optimum of ipopt
  IpoptWarmStart(const Ipopt::SmartPtr<TMINLP2TNLP> tnlp,
      SmartPtr<IpoptInteriorWarmStarter> warm_starter);
  /// Copy constructor
  IpoptWarmStart( const IpoptWarmStart &other, bool ownValues = 1);
  /// Abstract destructor
  virtual ~IpoptWarmStart();

  /// `Virtual constructor'
  virtual CoinWarmStart *clone() const
  {
    return new IpoptWarmStart(*this,1);
  }

  /** Generate the "differences" between two IpoptWarmStart.*/
  virtual CoinWarmStartDiff*
  generateDiff(const CoinWarmStart *const oldCWS) const;
  /** \brief Apply 'differences' to an Ipopt warm start.
   * What this actually does is get a copy to the vector of values stored
   in IpoptWarmStartDiff.*/
  virtual void
  applyDiff (const CoinWarmStartDiff *const cwsdDiff);
  /** Access to values_ vector. */
  const CoinPackedVector * values() const
  {
    if(tempValues_)
      return tempValues_;
    else
      return &values_;
  }
  /** Accessor to warm start information obecjt */
  SmartPtr<IpoptInteriorWarmStarter> warm_starter() const
  {
    return warm_starter_;
  }

  /// flush the starting point
  void flushPoint();

  ///Is this an empty warm start?
  bool empty() const
  {
    return empty_;
  }
private:
  /** Non zero values of the starting point. Primal and dual values are stored in the following order <p>
      <UL>
      <li> From 1 to CoinWarmStartBasis::numStrtucturals_ : values for primal variables (\f$ x \f$ ),
      <li> From numStructurals_+1 to  2numStructurals_ : values for dual variables associated to lower bound constraints on structurals (constraints \f$ l \leq x \f$).
      <li> From 2 numStructurals_+1 to  3 numStructurals_ : values for dual variables associated to lower bound constraints on structurals (constraints \f$ x \leq u\f$).
      <li> From 3 numStructurals_+1 to  3 numStructurals_ + numArtificials_ : values for dual varaibles associated with regular constraints (constraints \f$ g(x) = 0 \f$).
      </UL>
  */
  mutable CoinPackedVector values_;
  /** Temporary values not owned by this. */
  mutable CoinPackedVector * tempValues_;
  /** warm start information object */
  mutable SmartPtr<IpoptInteriorWarmStarter> warm_starter_;
  ///Say if warm start is empty
  bool empty_;
};

/** \brief Diff class for IpoptWarmStart.
 * Actually get the differences from CoinWarmStartBasis and stores the
 whole vector of values.
 \todo Find a way to free unused values.
*/
class IpoptWarmStartDiff : public CoinWarmStartBasisDiff
{
public:
  friend class IpoptWarmStart;
  /** Usefull constructor*/
  IpoptWarmStartDiff(CoinWarmStartBasisDiff * diff, const CoinPackedVector &values,
      SmartPtr<IpoptInteriorWarmStarter> warm_starter):
      CoinWarmStartBasisDiff(*diff),
      diffValues_(NULL),
      warm_starter_(NULL)//(warm_starter)
  {
    if(values.getNumElements()>0)
      diffValues_ = new CoinPackedVector(values);
  }
  /** Copy constructor. */
  IpoptWarmStartDiff(const IpoptWarmStartDiff &other):
      CoinWarmStartBasisDiff(other),
      diffValues_(NULL),
      warm_starter_(NULL)//other.warm_starter_)
  {
    if(other.diffValues_)
      diffValues_ = new CoinPackedVector(*other.diffValues_);
  }

  /// Abstract destructor
  virtual ~IpoptWarmStartDiff()
  {
    delete diffValues_;
  }

  /// `Virtual constructor'
  virtual CoinWarmStartDiff *clone() const
  {
    return new IpoptWarmStartDiff(*this);
  }

  /** Accessor to warm start information obecjt */
  SmartPtr<IpoptInteriorWarmStarter> warm_starter() const
  {
    return warm_starter_;
  }
  void flushPoint();
private:
  /** Values of the vector. */
  CoinPackedVector * diffValues_;

  /** warm start information object */
  SmartPtr<IpoptInteriorWarmStarter> warm_starter_;
};

}
#endif

#ifndef BonCouenneCbc_hpp
#define BonCouenneCbc_hpp

#include "BonCbcParam.hpp"
#include "BonCbc.hpp"



namespace Bonmin
{
  /** Class which performs optimization of an MINLP stored in an IpoptInterface. */
  class CouenneBab : public Bab
  {
  public:
    /** Constructor.*/
    CouenneBab();
    /** destructor.*/
    virtual ~CouenneBab();
    /** Perform a branch-and-bound on given IpoptInterface using passed parameters.*/
    virtual void branchAndBound(OsiTMINLPInterface * nlp,
        const BonminCbcParam&par);

    /**operator() performs the branchAndBound*/
    virtual void operator()(OsiTMINLPInterface * nlp, const BonminCbcParam& par)
    {
      CouenneBab::branchAndBound(nlp,par);
    }
  };
}
#endif

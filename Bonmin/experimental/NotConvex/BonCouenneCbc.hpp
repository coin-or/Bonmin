#ifndef BonminBB_hpp
#define BonminBB_hpp

#include "BonCbcParam.hpp"


class CbcObject;

namespace Bonmin
{
  class OsiTMINLPInterface;
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
      BonCouenneBab::branchAndBound(nlp,par);
    }

  };
}
#endif

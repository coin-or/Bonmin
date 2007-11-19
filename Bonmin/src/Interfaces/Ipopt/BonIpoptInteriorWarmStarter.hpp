// (C) Copyright International Business Machines Corporation 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter               IBM    2006-03-02

#ifndef __IPOPTINTERIORWARMSTARTER_HPP__
#define __IPOPTINTERIORWARMSTARTER_HPP__

#include "IpSmartPtr.hpp"
#include "IpNLP.hpp"
#include <vector>

using namespace Ipopt;
namespace Bonmin
{
  class IpoptInteriorWarmStarter : public ReferencedObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor. We give it the values of the current bounds so that
     *  it can figure out which variables are fixed for this NLP. */
    IpoptInteriorWarmStarter(Index n, const Number* x_L, const Number* x_u,
        Number nlp_lower_bound_inf,
        Number nlp_upper_bound_inf,
        bool store_several_iterates);

    /** Default destructor */
    ~IpoptInteriorWarmStarter();
    //@}

    /** Method for possibly storing another iterate during the current
     *  optimizatin for possible use for a warm start for a new
     *  problem */
    bool UpdateStoredIterates(AlgorithmMode mode,
        const IpoptData& ip_data,
        IpoptCalculatedQuantities& ip_cq);

    /** Method for doing whatever needs to be done after the parent NLP
     *  has been solved */
    bool Finalize();

    /** Method for computing the initial point based on the stored
     *  information */
    bool WarmStartIterate(Index n, const Number* x_l_new, const Number* x_u_new,
        IteratesVector& warm_start_iterate);

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default constructor. */
    IpoptInteriorWarmStarter();

    /** Copy Constructor */
    IpoptInteriorWarmStarter(const IpoptInteriorWarmStarter&);

    /** Overloaded Equals Operator */
    void operator=(const IpoptInteriorWarmStarter&);
    //@}

    //@{
    /** Value for a lower bound that denotes -infinity */
    Number nlp_lower_bound_inf_;
    /** Value for a upper bound that denotes infinity */
    Number nlp_upper_bound_inf_;
    /** Flag indicating whether more than one iterate is to be
     *  stored. */
    bool store_several_iterates_;
    //@}

    /** @name Copy of the bounds for the previously solved NLP.  This is
     *  required to find out the remapping for fixed variables, and it
     *  might also help to see how large the perturbation of the new
     *  problem is. */
    //@{
    Index n_;
    Number* x_l_prev_;
    Number* x_u_prev_;
    //@}

    /** @name Selected Iterates and quantities from the previous
     *  optimization */
    //@{
    Index n_stored_iterates_;
    std::vector<Index> stored_iter_;
    std::vector<SmartPtr<const IteratesVector> > stored_iterates_;
    std::vector<Number> stored_mu_;
    std::vector<Number> stored_nlp_error_;
    std::vector<Number> stored_primal_inf_;
    std::vector<Number> stored_dual_inf_;
    std::vector<Number> stored_compl_;
    //@}
  };
}
#endif

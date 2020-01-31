// (C) Copyright CNRS and others 2010
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, Université de la Méditérannée
// Hassan Hijazi, Orange Labs
//
// Date : 05/22/2010

#ifndef BonHeuristicInnerApproximation_HPP
#define BonHeuristicInnerApproximation_HPP
#include "BonOsiTMINLPInterface.hpp"
#include "BonBonminSetup.hpp"
#include "CbcHeuristic.hpp"
#include "CbcStrategy.hpp"

namespace Bonmin {
class SubMipSolver;
class HeuristicInnerApproximation: public CbcHeuristic {
public:

	/// Constructor with setup
	HeuristicInnerApproximation(BonminSetup * setup);

	/// Copy constructor
	HeuristicInnerApproximation(const HeuristicInnerApproximation &copy);

	/// Destructor
	~HeuristicInnerApproximation();

	/// Assignment operator
	HeuristicInnerApproximation & operator=(
			const HeuristicInnerApproximation & rhs);

	/// Clone
	virtual CbcHeuristic * clone() const {
		return new HeuristicInnerApproximation(*this);
	}

	/// Initialize method 
	void Initialize(BonminSetup * setup);

	/// Resets stuff if model changes
	virtual void resetModel(CbcModel * model) {
		setModel(model);
	}

	/// Performs heuristic
	virtual int solution(double &solutionValue, double *betterSolution);

	/** Register the options common to all local search based heuristics.*/
	static void registerOptions(
			Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions);

protected:
	/** Setup to use for local searches (will make copies).*/
	BonminSetup * setup_;

private:
	/// How often to do (code can change)
	int howOften_;

	/// A subsolver for MIP
	SubMipSolver * mip_;

        /// Number of Approximation points
        int nbAp_;

        void extractInnerApproximation(OsiTMINLPInterface & nlp, OsiSolverInterface &si,
                                       const double * x, bool getObj);

        bool getMyInnerApproximation(OsiTMINLPInterface &si, OsiCuts &cs, int ind,
                const double * x, const double * x2);


};
}

#endif

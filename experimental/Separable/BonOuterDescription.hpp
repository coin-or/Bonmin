// (C) Copyright CNRS and others 2010
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, Université de la Méditérannée
// Hassan Hijazi, Orange Labs
//
// Date : 05/22/2010

#ifndef OuterDescription_H
#define OuterDescription_H

#define INT_BIAS 0e-8

#include <string>
#include <iostream>

#include "BonOsiTMINLPInterface.hpp"
#include "CoinWarmStartBasis.hpp"

#include "BonCutStrengthener.hpp"

namespace Bonmin {

	/** Get the outer approximation constraints at provided point and only for the specified constraint (ind is the constraint or row number).
	 If x2 is different from NULL only add cuts violated by x2 by more than delta. **/
	void getMyOuterApproximation(OsiTMINLPInterface &si,
                                     OsiCuts &cs, int ind, const double * x,
			             int getObj, const double * x2,
                                     double theta, bool global);

/** Adds an outer description of problem to linear formulation.*/
void addOuterDescription(OsiTMINLPInterface &nlp, OsiSolverInterface &si,
		const double * x, int nbAp, bool getObj);

}
#endif

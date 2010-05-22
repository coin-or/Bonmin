// (C) Copyright CNRS and others 2010
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, Université de la Méditérannée
// Hassan Hijazi, Orange Labs
//
// Date : 05/22/2010

#include "BonOuterDescription.cpp"

/** Get the outer approximation constraints at provided point and only for the specified constraint 
 * (ind is the constraint or row number).
 * If x2 is different from NULL only add cuts violated by x2 by more than delta. **/
void OsiTMINLPInterface::getMyOuterApproximation(
                OsiTMINLPInterface &si, OsiCuts &cs, int ind,
		const double * x, int getObj, const double * x2, double theta,
		bool global) {
	int n, m, nnz_jac_g, nnz_h_lag;
	TNLP::IndexStyleEnum index_style;
	problem->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
	if (jRow_ == NULL || jCol_ == NULL || jValues_ == NULL)
		initializeJacobianArrays();
	assert(jRow_ != NULL);
	assert(jCol_ != NULL);
	vector<double> g(m);
	problem->eval_jac_g(n, x, 1, m, nnz_jac_g, NULL, NULL,
			jValues_);
	problem->eval_g(n, x, 1, m, g());
	//As jacobian is stored by cols fill OsiCuts with cuts
	CoinPackedVector cut;
	double lb;
	double ub;

	const double * rowLower = getRowLower();
	const double * rowUpper = getRowUpper();
	const double * colLower = getColLower();
	const double * colUpper = getColUpper();
	const double * duals = getRowPrice() + 2 * n;
	double infty = getInfinity();
	double nlp_infty = infty_;

	int rowIdx = ind;

	if (rowLower[rowIdx] > -nlp_infty)
		lb = rowLower[rowIdx] - g[rowIdx];
	else
		lb = -infty;
	if (rowUpper[rowIdx] < nlp_infty)
		ub = rowUpper[rowIdx] - g[rowIdx];
	else
		ub = infty;
	if (rowLower[rowIdx] > -infty && rowUpper[rowIdx] < infty) {
		if (duals[rowIdx] >= 0)// <= inequality
			lb = -infty;
		if (duals[rowIdx] <= 0)// >= inequality
			ub = infty;
	}

	for (int i = 0; i < nnz_jac_g; i++) {
		const int &rowIdx = jRow_[i];
		if (rowIdx == ind) {
			const int &colIdx = jCol_[i];
			//"clean" coefficient
			if (cleanNnz(jValues_[i], colLower[colIdx], colUpper[colIdx],
					rowLower[rowIdx], rowUpper[rowIdx], x[colIdx], lb, ub,
					tiny_, veryTiny_)) {
				cut.insert(colIdx, jValues_[i]);
				if (lb > -infty)
					lb += jValues_[i] * x[colIdx];
				if (ub < infty)
					ub += jValues_[i] * x[colIdx];
			}
		}
	}

	vector<int> cut2rowIdx(1);

	bool add = true;
	//Compute cut violation
	if (x2 != NULL) {
		double rhs = cut.dotProduct(x2);
		double violation = 0.;
		if (ub < infty)
			violation = std::max(violation, fabs(rhs - ub));
		if (lb > -infty)
			violation = std::max(violation, fabs(lb - rhs));
		if (violation < theta) {
			if (oaHandler_->logLevel() > 0)
				oaHandler_->message(CUT_NOT_VIOLATED_ENOUGH, oaMessages_)
						<< ind << violation << CoinMessageEol;
			add = false;
		}
		if (oaHandler_->logLevel() > 0)
			oaHandler_->message(VIOLATED_OA_CUT_GENERATED, oaMessages_) << ind
					<< violation << CoinMessageEol;
	} else if (oaHandler_->logLevel() > 0)
		oaHandler_->message(OA_CUT_GENERATED, oaMessages_) << ind
				<< CoinMessageEol;
	OsiRowCut newCut;
	//    if(lb[i]>-1e20) assert (ub[i]>1e20);

	if (IsValid(cutStrengthener_)) {
		const int& rowIdx = ind;
		bool retval = cutStrengthener_->ComputeCuts(cs, GetRawPtr(tminlp_),
				GetRawPtr(problem_), rowIdx, cut, lb, ub, g[rowIdx],
				rowLower[rowIdx], rowUpper[rowIdx], n, x, infty);
		if (!retval) {
			(*messageHandler()) << "error in cutStrengthener_->ComputeCuts\n";
			//exit(-2);
		}
	}
	if (add) {
		if (global) {
			newCut.setGloballyValidAsInteger(1);
		}
		//newCut.setEffectiveness(99.99e99);
		newCut.setLb(lb);
		newCut.setUb(ub);
		newCut.setRow(cut);
		if (oaHandler_->logLevel() > 2) {
			oaHandler_->print(newCut);
		}
		cs.insert(newCut);
	}
}

//END HASSAN




void addOuterDescription(OsiTMINLPInterface &nlp, OsiSolverInterface &si,
		const double * x, bool getObj) {
	int n;
	int m;
	int nnz_jac_g;
	int nnz_h_lag;
	TNLP::IndexStyleEnum index_style;
	//Get problem information
	problem->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);

	// Hassan OA initial description
	int OuterDesc = 0;
	app_->options()->GetEnumValue("initial_outer_description", OuterDesc,
			app_->prefix());
	if (OuterDesc == 0) {
		OsiCuts cs;

		double * p = CoinCopyOfArray(colLower, n);
		double * pp = CoinCopyOfArray(colLower, n);
		double * up = CoinCopyOfArray(colUpper, n);
#if 1
		int nbAp = 50;
		app_->options()->GetIntegerValue("number_approximations_initial_outer",
				nbAp, app_->prefix());
		std::vector<int> nbG(m, 0);// Number of generated points for each nonlinear constraint

		std::vector<double> step(n);

		for (int i = 0; i < n; i++) {

			if (colUpper[i] > 1e08) {
				up[i] = 0;
			}

			if (colUpper[i] > 1e08 || colLower[i] < -1e08 || (variableType[i]
					== TMINLP::BINARY) || (variableType[i] == TMINLP::INTEGER)) {
				step[i] = 0;
			} else
				step[i] = (up[i] - colLower[i]) / (nbAp * 20.);

			if (colLower[i] < -1e08) {
				p[i] = 0;
				pp[i] = 0;
			}
		}
		vector<double> g_p(m);
		vector<double> g_pp(m);
		for (int i = 0; (i < m); i++) {
                        if(constTypes_[i] != TNLP::NON_LINEAR) continue:
			getMyOuterApproximation(cs, i, p, 0, NULL, 10000, true);// Generate Tangents at current point    	 
		}
		for (int j = 1; j <= nbAp * 20; j++) {

			for (int i = 0; i < n; i++) {
				pp[i] += step[i];
			}

		}
		problem->eval_g(n, p, 1, m, g_p());
		problem->eval_g(n, pp, 1, m, g_pp());
		double diff = 0;
		int varInd = 0;
		for (int i = 0; (i < m); i++) {
                        if(constTypes_[i] != TNLP::NON_LINEAR) continue:
			if (varInd == n - 1)
				varInd = 0;
			diff = std::abs(g_p[i] - g_pp[i]);

			if (nbG[i] < nbAp && std::abs(g_p[i] - g_pp[i]) >= 0.005) {
				getMyOuterApproximation(cs, i, pp, 0, NULL, 10000, true);// Generate Tangents at current point
				p[varInd] = pp[varInd];
				nbG[i]++;
			}
			varInd++;
		}
		for (int i = 0; i < m ; i++) {
                        if(constTypes_[i] != TNLP::NON_LINEAR) continue:
			getMyOuterApproximation(cs, i, up, 0, NULL, 10000, true);// Generate Tangents at current point
		}
		si.applyCuts(cs);
		delete [] p;
		delete [] pp;
		delete [] up;
	}

}


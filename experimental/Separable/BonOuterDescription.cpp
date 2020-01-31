// (C) Copyright CNRS and others 2010
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Pierre Bonami, Université de la Méditérannée
// Hassan Hijazi, Orange Labs
//
// Date : 05/22/2010

#include "BonOuterDescription.hpp"
#include "BonOsiTMINLPInterface.hpp"

namespace Bonmin{


//Copied from OsiTMINLPInterface

//A procedure to try to remove small coefficients in OA cuts (or make it non small
static inline
bool cleanNnz(double &value, double colLower, double colUpper,
    double rowLower, double rowUpper, double colsol,
    double & lb, double &ub, double tiny, double veryTiny)
{
  if(fabs(value)>= tiny) return 1;

  if(fabs(value)<veryTiny) return 0;//Take the risk?

  //try and remove
  double infty = 1e20;
  bool colUpBounded = colUpper < 10000;
  bool colLoBounded = colLower > -10000;
  bool rowNotLoBounded =  rowLower <= - infty;
  bool rowNotUpBounded = rowUpper >= infty;
  bool pos =  value > 0;

  if(colLoBounded && pos && rowNotUpBounded) {
    lb += value * (colsol - colLower);
    return 0;
  }
  else
    if(colLoBounded && !pos && rowNotLoBounded) {
      ub += value * (colsol - colLower);
      return 0;
    }
    else
      if(colUpBounded && !pos && rowNotUpBounded) {
        lb += value * (colsol - colUpper);
        return 0;
      }
      else
        if(colUpBounded && pos && rowNotLoBounded) {
          ub += value * (colsol - colUpper);
          return 0;
        }
  //can not remove coefficient increase it to smallest non zero
  if(pos) value = tiny;
  else
    value = - tiny;
  return 1;
}


/** Get the outer approximation constraints at provided point and only for the specified constraint 
 * (ind is the constraint or row number).
 * If x2 is different from NULL only add cuts violated by x2 by more than delta. **/
void getMyOuterApproximation(
                OsiTMINLPInterface &si, OsiCuts &cs, int ind,
		const double * x, int getObj, const double * x2, double theta,
		bool global) {
	int n, m, nnz_jac_g, nnz_h_lag;
	Ipopt::TNLP::IndexStyleEnum index_style;
        TMINLP2TNLP* problem = si.problem();
	problem->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);

        double g_i = 0;
	problem->eval_gi(n, x, 1, ind, g_i);
        vector<int> jCol(n);
        int nnz;
        problem->eval_grad_gi(n, x, 0, ind, nnz, jCol(), NULL);
        vector<double> jValues(nnz);
        problem->eval_grad_gi(n, x, 0, ind, nnz, NULL, jValues());
	//As jacobian is stored by cols fill OsiCuts with cuts
	CoinPackedVector cut;
	double lb;
	double ub;

	const double * rowLower = si.getRowLower();
	const double * rowUpper = si.getRowUpper();
	const double * colLower = si.getColLower();
	const double * colUpper = si.getColUpper();
	const double * duals = si.getRowPrice() + 2 * n;
	double infty = si.getInfinity();
	double nlp_infty = infty;
	int rowIdx = ind;

	if (rowLower[rowIdx] > -nlp_infty)
		lb = rowLower[rowIdx] - g_i;
	else
		lb = -infty;
	if (rowUpper[rowIdx] < nlp_infty)
		ub = rowUpper[rowIdx] - g_i;
	else
		ub = infty;
	if (rowLower[rowIdx] > -infty && rowUpper[rowIdx] < infty) {
		if (duals[rowIdx] >= 0)// <= inequality
			lb = -infty;
		if (duals[rowIdx] <= 0)// >= inequality
			ub = infty;
	}

        double tiny = 1e-08;
        double veryTiny = 1e-20;

	for (int i = 0; i < nnz; i++) {
	  if(index_style == Ipopt::TNLP::FORTRAN_STYLE) jCol[i]--;
	  const int &colIdx = jCol[i];
	//"clean" coefficient
	if (cleanNnz(jValues[i], colLower[colIdx], colUpper[colIdx],
			rowLower[rowIdx], rowUpper[rowIdx], x[colIdx], lb, ub,
			tiny, veryTiny)) {
		cut.insert(colIdx, jValues[i]);
        
		if (lb > -infty)
		lb += jValues[i] * x[colIdx];
		if (ub < infty)
		ub += jValues[i] * x[colIdx];
	}
	}


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
			add = false;
		}
	} 
	OsiRowCut newCut;
	//    if(lb[i]>-1e20) assert (ub[i]>1e20);

	if (add) {
		if (global) {
		newCut.setGloballyValidAsInteger(1);
		}
		//newCut.setEffectiveness(99.99e99);

        //********* Perspective Extension ********//
        const int* ids = problem->get_const_xtra_id(); // vector of indices corresponding to the binary variable activating the corresponding constraint
        // Get the index of the corresponding indicator binary variable
        int binary_id = (ids == NULL) ? -1 : ids[ind];// index corresponding to the binary variable activating the corresponding constraint
        if(binary_id>0) {// If this hyperplane is a linearization of a disjunctive constraint, we link its righthand (or lefthand) side to the corresponding indicator binary variable
            if (lb > -infty) { // ∂x ≥ lb => ∂x - lb*z ≥ 0
                cut.insert(binary_id, -lb);
                newCut.setLb(0);
                newCut.setUb(ub);
                
            }
            if (ub < infty) { // ∂x ≤ ub => ∂x - ub*z ≤ 0
                cut.insert(binary_id, -ub);
                newCut.setLb(lb);
                newCut.setUb(0);
                
            }
        }
        else {
            newCut.setLb(lb);
            newCut.setUb(ub);
        }
        //********* Perspective Extension ********//

		newCut.setRow(cut);
		//newCut.print();
		cs.insert(newCut);
	}
}
    
/* Old Uniform sampling method
void addOuterDescription(OsiTMINLPInterface &nlp, OsiSolverInterface &si,
		const double * x, int nbAp, bool getObj) {
	int n;
	int m;
	int nnz_jac_g;
	int nnz_h_lag;
	Ipopt::TNLP::IndexStyleEnum index_style;
	//Get problem information
        TMINLP2TNLP* problem = nlp.problem();
	problem->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);

   const double * colLower = nlp.getColLower();
   const double * colUpper = nlp.getColUpper();
   const Bonmin::TMINLP::VariableType* variableType = problem->var_types();
   vector<Ipopt::TNLP::LinearityType>  constTypes(m);
   problem->get_constraints_linearity(m, constTypes());
	// Hassan OA initial description
	int OuterDesc = 0;
	if (OuterDesc == 0) {
		OsiCuts cs;

		double * p = CoinCopyOfArray(nlp.getColLower(), n);
		double * pp = CoinCopyOfArray(nlp.getColLower(), n);
		double * up = CoinCopyOfArray(nlp.getColUpper(), n);
		//b->options()->GetIntegerValue("number_approximations_initial_outer",nbAp, b->prefix());
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
				step[i] = (up[i] - colLower[i]) / (nbAp);

			if (colLower[i] < -1e08) {
				p[i] = 0;
				pp[i] = 0;
			}
		}
		vector<double> g_p(m);
		vector<double> g_pp(m);
		for (int i = 0; (i < m); i++) {
                        if(constTypes[i] != Ipopt::TNLP::NON_LINEAR) continue;
			getMyOuterApproximation(nlp, cs, i, p, 0, NULL, 10000, true);// Generate Tangents at current point    	 
		}
		for (int j = 1; j <= nbAp; j++) {

			for (int i = 0; i < n; i++) {
				pp[i] += step[i];
			}

		
		problem->eval_g(n, p, 1, m, g_p());
		problem->eval_g(n, pp, 1, m, g_pp());
		double diff = 0;
		int varInd = 0;
		for (int i = 0; (i < m); i++) {
                        if(constTypes[i] != Ipopt::TNLP::NON_LINEAR) continue;
			if (varInd == n - 1)
				varInd = 0;
			diff = std::abs(g_p[i] - g_pp[i]);

			if (nbG[i] < nbAp && diff ) {
				getMyOuterApproximation(nlp, cs, i, pp, 0, NULL, 10000, true);// Generate Tangents at current point
				p[varInd] = pp[varInd];
				nbG[i]++;
			}
			varInd++;
		}
		}
		for (int i = 0; i < m ; i++) {
                        if(constTypes[i] != Ipopt::TNLP::NON_LINEAR) continue;
			getMyOuterApproximation(nlp, cs, i, up, 0, NULL, 10000, true);// Generate Tangents at current point
		}
		
		si.applyCuts(cs);
		delete [] p;
		delete [] pp;
		delete [] up;
	}

} */
    // New curvature based sampling method
    void addOuterDescription(OsiTMINLPInterface &nlp, OsiSolverInterface &si,
                             const double * x, int nbAp, bool getObj) {
        int n;
        int m;
        int nnz_jac_g;
        int nnz_h_lag;
        Ipopt::TNLP::IndexStyleEnum index_style;
        //Get problem information
        TMINLP2TNLP* problem = nlp.problem();
        problem->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
        
        const double * colLower = nlp.getColLower();
        const double * colUpper = nlp.getColUpper();
        const Bonmin::TMINLP::VariableType* variableType = problem->var_types();
        vector<Ipopt::TNLP::LinearityType>  constTypes(m);
        problem->get_constraints_linearity(m, constTypes());
        // Hassan OA initial description
        int OuterDesc = 0;
        if (OuterDesc == 0) {
            OsiCuts cs;
            
            double * p = CoinCopyOfArray(nlp.getColLower(), n);
            double * pp = CoinCopyOfArray(nlp.getColLower(), n);
            double * up = CoinCopyOfArray(nlp.getColUpper(), n);
            //b->options()->GetIntegerValue("number_approximations_initial_outer",nbAp, b->prefix());
            std::vector<int> nbG(m, 2);// Number of generated points for each nonlinear constraint, we always generate two cuts at lower and upper bounds.
            
            std::vector<double> step(n);
            
            for (int i = 0; i < n; i++) {
                
                if (colUpper[i] > 1e08) {
                    up[i] = 0;
                }
                
                if (colUpper[i] > 1e08 || colLower[i] < -1e08 || (variableType[i]
                                                                  == TMINLP::BINARY) || (variableType[i] == TMINLP::INTEGER)) {
                    step[i] = 0;
                } else
                    step[i] = (up[i] - colLower[i]) / 2e02; // if step[i]!=0 then variable i appears in one univariate nonlinear function, a small step is computed in order to perform curavture based sampling
                
                if (colLower[i] < -1e08) {
                    p[i] = 0;
                    pp[i] = 0;
                }
            }
            vector<double> g_p(m);
            double g_p_i, g_pp_i;
            problem->eval_g(n, p, 1, m, g_p()); // Evaluate function g at lowerbounds
            vector<double> g_pp(m);
            vector<double> g_up(m);
            problem->eval_g(n, up, 1, m, g_up()); // Evaluate function g at upperbounds
            
            for (int i = 0; (i < m); i++) {
                if(constTypes[i] != Ipopt::TNLP::NON_LINEAR) continue;
                getMyOuterApproximation(nlp, cs, i, p, 0, NULL, 10000, true);// Generate Tangents at lowerbounds    	 
            }
            vector<double> thr(m); // minimum threshold value for the variation of function g (curvature) for generating a cut at point pp
            for (int i = 0; i < m; i++) {
                thr[i] = std::abs(g_up[i]-g_p[i])/nbAp;
            }
            double diff = 0;
            for (int i = 0; (i < m); i++) { // Generate Outer-Approximation initial cuts for all nonlinear constraints
                if(constTypes[i] != Ipopt::TNLP::NON_LINEAR) continue;
                p = CoinCopyOfArray(nlp.getColLower(), n); // For each nonlinear constraint reset all points to lowerbounds
                pp = CoinCopyOfArray(nlp.getColLower(), n);
                while (nbG[i] < nbAp) { // Iterate untill the number of initial approximations is reached for each nonlinear constraint
                
                    // Curvature sampling
                    for (int j = 0; j < n; j++) {
                        pp[j] += step[j];
                    }
                    problem->eval_gi(n, p, 1, i, g_p_i);
                    problem->eval_gi(n, pp, 1, i, g_pp_i);
                    diff = std::abs(g_p_i - g_pp_i);  
                        
                    if (diff>=thr[i] ) {
                        getMyOuterApproximation(nlp, cs, i, pp, 0, NULL, 10000, true);// Generate Tangents at current point
                        for (int j = 0; j < n; j++) {
                            p[j] = pp[j]; // Move all previous points to the current one
                        }
                        
                        nbG[i]++;
                    }
                }
                
            }
            
            for (int i = 0; i < m ; i++) {
                if(constTypes[i] != Ipopt::TNLP::NON_LINEAR) continue;
                getMyOuterApproximation(nlp, cs, i, up, 0, NULL, 10000, true);// Generate Tangents at upperbounds
            }
            
            si.applyCuts(cs);
            delete [] p;
            delete [] pp;
            delete [] up;
        
        }


    }
}


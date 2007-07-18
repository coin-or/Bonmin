/*
 * Name:    alphaConvexify.cpp
 * Authors: Pierre Bonami
 *          Stefan Vigerske
 *          Pietro Belotti
 * Purpose: create alpha-convexification of a quadratic expression
 *
 * (C) Pietro Belotti 2007. This file is licensed under the Common Public License (CPL)
 */

#include <exprQuad.hpp>

#include <OsiSolverInterface.hpp>
#include <IpLapack.hpp>

  /** Wrapper for LAPACK subroutine DSYEV.  Compute the Eigenvalue
   *  decomposition for a given matrix.  If compute_eigenvectors is
   *  true, a will contain the eigenvectors in its columns on
   *  return.  */
//  void IpLapackDsyev(bool compute_eigenvectors, Index ndim, Number *a,
//                     Index lda, Number *w, Index& info);


/// [Stefan] fills in dCoeffLo_, dCoeffUp_, and dIndex_ for the convex
/// under- and overestimator of this expression
/** Computes alpha coefficients for an alpha under- and overestimator of the quadratic term.
 * For the underestimator, dCoeffLo_ is computed such that
 * x^TQx + sum_i dCoeffLo_i (x_i - lb x_i)(ub x_i - x-i) is convex and underestimating (alpha_i is negative),
 * Regarding the overestimator, dCoeffUp_ are computed such that
 * x^TQx + sum_i dCoeffUp_i (x_i - lb x_i)(ub x_i - x-i) is concave and overestimating (alpha_i is positive).
 * If the method hasn't been called, dIndex_ will be NULL.
 * If Q is positive semidefinite, then dCoeffLo_ will be NULL.
 * If Q is negative semidefinite, then dCoeffUp_ will be NULL.
*/
void exprQuad::alphaConvexify (const OsiSolverInterface &si) {
	if (getnQTerms()==0) {
		nDiag_=0;
		return;
	}

	// inverse of dIndex_ mapping, for each variable tell me the index that it will have in dIndex_, or -1 if not there
	int* indexmap=new int[si.getNumCols()];
	for (int i=0; i<si.getNumCols(); ++i)
		indexmap[i]=-1;
	if (dIndex_==NULL) { // first time called... check which variables are there, and where we will put in the dIndex_ array
		nDiag_=0;
		for (int i=0; i<getnQTerms(); ++i) {
			if (indexmap[getQIndexI()[i]]==-1) {
				indexmap[getQIndexI()[i]]=nDiag_;
				++nDiag_;
			}
			if (indexmap[getQIndexJ()[i]]==-1) {
				indexmap[getQIndexJ()[i]]=nDiag_;
				++nDiag_;
			}
		}
		dIndex_=new int[nDiag_];
		for (int i=0; i<si.getNumCols(); ++i) {
			if (indexmap[i]>-1) {
				dIndex_[indexmap[i]]=i;
			}
		}
	} else {
		for (int i=0; i<nDiag_; ++i)
			indexmap[dIndex_[i]]=i;
	}

	// box diameter
	double* diam=new double[nDiag_];
	for (int i=0; i<nDiag_; ++i)
		diam[i]=si.getColUpper()[dIndex_[i]]-si.getColLower()[dIndex_[i]];

	// lower triangular of quadratic term matrix, scaled by box diameter
	double* matrix=new double[nDiag_*(nDiag_+1)/2];
	for (int i=0; i<(nDiag_*(nDiag_+1))/2; ++i)
		matrix[i]=0.;
	for (int i=0; i<getnQTerms(); ++i) {
		int row=indexmap[getQIndexI()[i]];
		int col=indexmap[getQIndexJ()[i]];
		matrix[((row+1)*row)/2+col]=getQCoeffsI()[i]*diam[row]*diam[col];
	}

	// compute minimum and maximum eigenvalue of matrix
	// ok, currently computes all eigenvalues
	double* eigval=new double[nDiag_];
	int info;
	Ipopt::IpLapackDsyev(false, nDiag_, matrix, nDiag_, eigval, info);
	if (info!=0) {
		//TODO error handling
	}

	// if min. eigenvalue negative, setup dCoeffLo_
	if (eigval[0]<0) {
		if (dCoeffLo_==NULL)
			dCoeffLo_ = new CouNumber[nDiag_];
		for (int i=0; i<nDiag_; ++i) {
			if (diam[i]==0.)
				dCoeffLo_[i] = 0.;
			else
				dCoeffLo_[i] = eigval[0]/(diam[i]*diam[i]);
		}
	} else { // quadratic term is convex, no convexification needed
		if (dCoeffLo_)
			delete dCoeffLo_;
	}

	// if max. eigenvalue is positive, setup dCoeffUp_
	if (eigval[nDiag_-1]>0) {
		if (dCoeffUp_==NULL)
			dCoeffUp_ = new CouNumber[nDiag_];
		for (int i=0; i<nDiag_; ++i) {
			if (diam[i]==0.)
				dCoeffUp_[i] = 0.;
			else
				dCoeffUp_[i] = eigval[nDiag_-1]/(diam[i]*diam[i]);
		}
	} else { // quadratic term is concave, no "concavification" needed
		if (dCoeffUp_)
			delete dCoeffUp_;
	}

	delete[] matrix;
	delete[] diam;
	delete[] eigval;
}

/*
 * Name:    genSDPcut.cpp
 * Author:  Pietro Belotti
 * Purpose: generate cut of the form v^T X' w >= 0, where v and w 
 *          form the coefficient set and
 *          
 *         / 1  x' \
 *  X' =   \ x  X  /
 *
 * is the set of variables
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License.
 */

#include <OsiCuts.hpp>
#include <OsiRowCut.hpp>
#include <SdpCutGen.hpp>


/// use v1 and v2 to create cut v1^T X' v2 >= 0

void SdpCutGen::genSDPcut (OsiCuts &cs, double *v1, double *v2) const {

  int N      = n_*(n_+3)/2, // total number of terms in the cut: n(n+1)/2 for X and n for x
      nterms = 0, 
      np     = n_+1;

  OsiRowCut *cut   = new OsiRowCut;
  double    *coeff = new double [N];
  int       *ind   = new int    [N];

  /// coefficients for X_ij

  for (int i=1; i<np; i++)
    for (int j=i; j<np; j++) {

      double coeff0 = 
	v1 [i] * v2 [j] + 
	v1 [j] * v2 [i];

      if (fabs (coeff0) > 1e-6) {
	coeff [nterms] = (i==j) ? (0.5 * coeff0) : (coeff0);
	ind   [nterms++] = indexQ (i-1, j-1, n_);
      }
    }

  /// coefficients for x_i

  for (int i=1; i<np; i++) {

    double coeff0 = 
      v1 [i] * v2 [0] + 
      v1 [0] * v2 [i];

    if (fabs (coeff0) > 1e-6) {
      coeff [nterms]   = coeff0;
      ind   [nterms++] = i-1;
    }
  }

  cut -> setRow (nterms, ind, coeff);
  cut -> setLb (- *v1 * *v2);

  cs.insert (cut);

  delete cut;
  delete [] ind;
  delete [] coeff;
}

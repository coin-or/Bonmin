#include <OsiCuts.hpp>
#include <OsiRowCut.hpp>
#include <SdpCutGen.hpp>

void SdpCutGen::genSDPcut (OsiCuts &cs, double *x1, double *x2) const {

  int N      = n_*(n_+3)/2, 
      nterms = 0, 
      np     = n_+1;

  OsiRowCut *cut   = new OsiRowCut;
  double    *coeff = new double [N];
  int       *ind   = new int    [N];

  for (int i=1; i<np; i++)
    for (int j=i; j<np; j++) {

      double coeff0 = 
	x1 [i] * x2 [j] + 
	x1 [j] * x2 [i];

      if (fabs (coeff0) > 1e-6) {
	coeff [nterms] = (i==j) ? (0.5 * coeff0) : (coeff0);
	ind   [nterms++] = indexQ (i-1, j-1, n_);
      }
    }

  for (int i=1; i<np; i++) {

    double coeff0 = 
      x1 [i] * x2 [0] + 
      x1 [0] * x2 [i];

    if (fabs (coeff0) > 1e-6) {
      coeff [nterms]   = coeff0;
      ind   [nterms++] = i-1;
    }
  }

  cut -> setRow (nterms, ind, coeff);
  cut -> setLb (- *x1 * *x2);

  cs.insert (cut);

  delete cut;
  delete [] ind;
  delete [] coeff;
}
